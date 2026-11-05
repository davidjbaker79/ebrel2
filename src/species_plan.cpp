//-------------------- Pre-compute species plan of modelling -------------------

#include <cmath>
#include <algorithm>
#include <vector>

#include "species_plan.h"

//-------------------------- Local helper functions ----------------------------

namespace {

// Get number of cells
inline std::size_t n_cells(int dim_x, int dim_y) {
  return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
}

}

//----------------------------- Main functions ---------------------------------

// ---- Function to build species' plan (e.g., seed, roi etc) ----
SpeciesPlan build_species_plan(
    const std::vector<double>& SD,
    const std::vector<double>& SxH,   // R column-major: h*n_species + sp
    const std::vector<int>&    D,     // disp_raw per species
    const std::vector<double>& O,     // targets
    const std::vector<float>& action_weight,  // length n_actions
    double marginal_occupancy_weight,
    int n_actions, int n_species, int n_habitats,
    int dim_x, int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
    const std::vector<uint8_t>& LM,
    const std::vector<int8_t>& E,
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c
) {
  const std::size_t cells = n_cells(dim_x, dim_y);

  if (SD.size() != static_cast<std::size_t>(n_species) * cells)
    throw std::runtime_error("SD size mismatch in build_species_plan_soa");
  if (SxH.size() != static_cast<std::size_t>(n_species) * static_cast<std::size_t>(n_habitats))
    throw std::runtime_error("SxH size mismatch in build_species_plan_soa");
  if (D.size() != static_cast<std::size_t>(n_species))
    throw std::runtime_error("D size mismatch in build_species_plan_soa");
  if (action_weight.size() != static_cast<std::size_t>(n_actions))
    throw std::runtime_error("Action_weight must have length n_actions");
  if (O.size() != static_cast<std::size_t>(n_species))
    throw std::runtime_error("O size mismatch in build_species_plan");
  if (n_actions != n_habitats + 1)
    throw std::runtime_error(
        "Expected n_actions == n_habitats + 1 "
        "(final action is improvement)"
    );
  if (LM.size() != cells)
    throw std::runtime_error(
        "LM size mismatch in build_species_plan"
    );
  
  if (!std::isfinite(marginal_occupancy_weight) ||
      marginal_occupancy_weight <= 0.0 ||
      marginal_occupancy_weight > 1.0) {
    throw std::runtime_error(
        "marginal_occupancy_weight must be finite and in (0, 1]"
    );
  }
  
  if (E.size() != cells)
    throw std::runtime_error(
        "E size mismatch in build_species_plan"
    );
  
  SpeciesPlan P;
  P.n_species = n_species;
  P.n_actions = n_actions;
  P.n_habitats = n_habitats;
  
  P.active.assign(n_species, 0u);
  P.fast_path.assign(n_species, 0u);
  P.disp.assign(n_species, 0);
  P.r0.assign(n_species, 0); P.r1.assign(n_species, -1);
  P.c0.assign(n_species, 0); P.c1.assign(n_species, -1);
  P.W.assign(n_species, 0);  P.H.assign(n_species, 0);
  P.O_n.assign(n_species, 0.0);

  // ---- Action weights (global, length n_action + 1)
  P.action_weight = action_weight;
  
  // ---- Marginal occupancy weights
  P.marginal_occupancy_weight = marginal_occupancy_weight;
  
  // ---- Ns: occupied cells per species from SD
  P.Ns.assign(n_species, 0);
  for (int sp = 0; sp < n_species; ++sp) {
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;
    int count = 0;
    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] > 0.0 && LM[k])
        count++;
    }
    P.Ns[sp] = count;
  }
  
  P.suitable.assign(static_cast<std::size_t>(n_species) * static_cast<std::size_t>(n_habitats), 0u);

  P.seed_start.assign(n_species, 0u);
  P.seed_len.assign(n_species, 0u);
  
  P.seed_weight_sum.assign(n_species, 0.0);
  
  P.seed_pool.clear();
  P.seed_weight.clear();
  
  const std::size_t reserve_n =
    1024u * static_cast<std::size_t>(n_species);
  
  P.seed_pool.reserve(reserve_n);
  P.seed_weight.reserve(reserve_n);
  
  for (int sp = 0; sp < n_species; ++sp) {

    // ---- Suitable flags for this species (contiguous in P.suitable)
    const std::size_t base_h =
      static_cast<std::size_t>(sp) *
      static_cast<std::size_t>(n_habitats);
    
    for (int h = 0; h < n_habitats; ++h) {
      
      // Incoming SxH is an R species × habitat matrix flattened column-major:
      //   index = habitat * n_species + species
      const std::size_t sxh_idx =
        static_cast<std::size_t>(h) *
        static_cast<std::size_t>(n_species) +
        static_cast<std::size_t>(sp);
      
      // Store internally species-major for fast downstream access:
      //   index = species * n_habitats + habitat
      P.suitable[
      base_h + static_cast<std::size_t>(h)
      ] =
        (SxH[sxh_idx] > 0.0) ? 1u : 0u;
    }

    const int disp_raw = D[sp];
    const int disp = std::min(disp_raw, universal_disp_thres);
    P.disp[sp] = static_cast<int16_t>(disp);
    P.fast_path[sp] = (disp_raw >= universal_disp_thres) ? 1u : 0u;

    // ---- Seeds + bbox (but store seeds into global flat pool)
    // Any observed occurrence is treated as evidence that suitable habitat
    // exists somewhere within the planning unit. Therefore seed status does
    // not depend on compatibility between the coarse existing-habitat class E
    // and the species-habitat association SxH.
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;

    int min_r = dim_y, max_r = -1;
    int min_c = dim_x, max_c = -1;

    const uint32_t start = static_cast<uint32_t>(P.seed_pool.size());
    P.seed_start[sp] = start;

    for (std::size_t k = 0; k < cells; ++k) {
      
      // Must be an observed occurrence on a valid landscape cell.
      if (SD[base_sd + k] <= 0.0) continue;
      if (!LM[k]) continue;
      
      // Existing mapped habitat at this cell.
      const int he = E[k];
      
      // Does the mapped main habitat match one of the species'
      // suitable habitats?
      const bool mapped_habitat_suitable =
        he >= 0 &&
        he < n_habitats &&
        P.suitable[
      base_h + static_cast<std::size_t>(he)
        ] != 0u;
      
      // Full cell-equivalent if occurrence agrees with mapped habitat.
      // Otherwise count it as marginal/sub-grid occupancy.
      const float seed_w =
        mapped_habitat_suitable
        ? 1.0f
      : static_cast<float>(marginal_occupancy_weight);
      
      // Store seed and its corresponding occupancy weight.
      P.seed_pool.push_back(
        static_cast<SeedT>(k)
      );
      
      P.seed_weight.push_back(seed_w);
      
      P.seed_weight_sum[sp] +=
        static_cast<double>(seed_w);
      
      // Bounding box still uses every observed occurrence,
      // irrespective of its occupancy weight.
      const int r = cell_r[k];
      const int c = cell_c[k];
      
      if (r < min_r) min_r = r;
      if (r > max_r) max_r = r;
      if (c < min_c) min_c = c;
      if (c > max_c) max_c = c;
    }

    // ---- Number of seeds
    const uint32_t len = static_cast<uint32_t>(P.seed_pool.size()) - start;
    P.seed_len[sp] = len;

    // ---- Absolute target in weighted cell-equivalent units.
    //
    // Baseline distribution size is the sum of seed weights:
    //   1.0 = occurrence in suitable mapped habitat
    //   marginal_occupancy_weight = occurrence in unsuitable mapped habitat
    //
    // O[sp] is the required proportional gain relative to this
    // weighted baseline distribution.
    
    const double o = O[static_cast<std::size_t>(sp)];
    
    if (!std::isfinite(o) || o < 0.0) {
      throw std::runtime_error(
          "Target O must be finite and non-negative in build_species_plan."
      );
    }
    
    P.O_n[sp] = P.seed_weight_sum[sp] * o;

    // ---- Mark inactive species
    if (len == 0u) {
      // inactive species
      P.active[sp] = 0u;
      continue;
    }

    // ---- ROI expand
    const int R_raw = (max_disp_steps > 0)
      ? (disp_raw * max_disp_steps)
      : std::max(dim_x, dim_y);
    const int R = std::min(R_raw, roi_cap);

    int r0 = std::max(0,         min_r - R);
    int r1 = std::min(dim_y - 1, max_r + R);
    int c0 = std::max(0,         min_c - R);
    int c1 = std::min(dim_x - 1, max_c + R);

    // tighten vertically using row spans
    while (r0 <= r1 && (row_last_land[r0] < c0 || row_first_land[r0] > c1)) ++r0;
    while (r1 >= r0 && (row_last_land[r1] < c0 || row_first_land[r1] > c1)) --r1;
    if (r0 > r1) { P.active[sp] = 0u; continue; }

    // tighten horizontally using col spans
    while (c0 <= c1 && (col_last_land[c0] < r0 || col_first_land[c0] > r1)) ++c0;
    while (c1 >= c0 && (col_last_land[c1] < r0 || col_first_land[c1] > r1)) --c1;
    if (c0 > c1) { P.active[sp] = 0u; continue; }

    const int W = c1 - c0 + 1;
    const int H = r1 - r0 + 1;
    if (W <= 0 || H <= 0) { P.active[sp] = 0u; continue; }

    P.r0[sp] = r0; P.r1[sp] = r1;
    P.c0[sp] = c0; P.c1[sp] = c1;
    P.W[sp]  = W;  P.H[sp]  = H;
    P.active[sp] = 1u;
  }
  
  // ---- Reverse index: for each cell, which species have that cell as a seed
  P.cell_species_len.assign(cells, 0u);
  
  // Count memberships
  for (int sp = 0; sp < n_species; ++sp) {
    const uint32_t start = P.seed_start[sp];
    const uint32_t len   = P.seed_len[sp];
    
    for (uint32_t i = 0; i < len; ++i) {
      const uint32_t cell = static_cast<uint32_t>(P.seed_pool[start + i]);
      P.cell_species_len[cell]++;
    }
  }
  
  // Prefix sums -> starts
  P.cell_species_start.resize(cells, 0u);
  uint32_t total_memberships = 0u;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    P.cell_species_start[cell] = total_memberships;
    total_memberships += P.cell_species_len[cell];
  }
  
  // ---- Fill flattened reverse pools
  P.cell_species_pool.assign(
    total_memberships,
    -1
  );
  
  std::vector<uint32_t> cursor =
    P.cell_species_start;
  
  for (int sp = 0; sp < n_species; ++sp) {
    
    const uint32_t start =
      P.seed_start[sp];
    
    const uint32_t len =
      P.seed_len[sp];
    
    for (uint32_t i = 0; i < len; ++i) {
      
      const uint32_t seed_idx =
        start + i;
      
      const uint32_t cell =
        static_cast<uint32_t>(
          P.seed_pool[seed_idx]
        );
      
      const uint32_t pos =
        cursor[cell]++;
      
      P.cell_species_pool[pos] =
        sp;
      
    }
  }
  
  if (P.seed_weight.size() != P.seed_pool.size()) {
    throw std::runtime_error(
        "seed_weight and seed_pool size mismatch in build_species_plan"
    );
  }

  return P;
}
