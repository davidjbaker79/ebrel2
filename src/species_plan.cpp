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
    const std::vector<double>& SxH,   // row-major sp*n_habitats + h
    const std::vector<int>&    D,     // disp_raw per species
    const std::vector<double>& O,     // targets
    const std::vector<float>& action_weight,  // length n_actions
    int n_actions, int n_species, int n_habitats,
    int dim_x, int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
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
  
  // ---- Ns: occupied cells per species from SD
  P.Ns.assign(n_species, 0);
  for (int sp = 0; sp < n_species; ++sp) {
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;
    int count = 0;
    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] > 0.0) count++;
    }
    P.Ns[sp] = count;
  }
  
  P.suitable.assign(static_cast<std::size_t>(n_species) * static_cast<std::size_t>(n_habitats), 0u);

  P.seed_start.assign(n_species, 0u);
  P.seed_len.assign(n_species, 0u);
  P.seed_pool.clear();
  P.seed_pool.reserve(1024u * static_cast<std::size_t>(n_species)); // heuristic
  
  for (int sp = 0; sp < n_species; ++sp) {

    // ---- Suitable flags for this species (contiguous in P.suitable)
    const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_habitats);
    for (int h = 0; h < n_habitats; ++h) {
      P.suitable[base_h + static_cast<std::size_t>(h)] =
        (SxH[base_h + static_cast<std::size_t>(h)] > 0.0) ? 1u : 0u;
    }

    const int disp_raw = D[sp];
    const int disp = std::min(disp_raw, universal_disp_thres);
    P.disp[sp] = static_cast<int16_t>(disp);
    P.fast_path[sp] = (disp_raw >= universal_disp_thres) ? 1u : 0u;

    // ---- Seeds + bbox (but store seeds into global flat pool)
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;

    int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;

    const uint32_t start = static_cast<uint32_t>(P.seed_pool.size());
    P.seed_start[sp] = start;

    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] != 1.0) continue;

      const int ehx = E[k];
      if (ehx == -1) continue;
      if (!P.suitable[base_h + static_cast<std::size_t>(ehx)]) continue;
      
      P.seed_pool.push_back(static_cast<SeedT>(k));

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

    // ---- Absolute target in cell-equivalent contribution units.
    //
    // O[sp] is the required proportional gain relative to the species'
    // baseline distribution. The target is deliberately allowed to be
    // fractional because action_weight is a per-cell contribution and may
    // itself be fractional.
    const double o = O[static_cast<std::size_t>(sp)];
    
    if (!std::isfinite(o) || o < 0.0) {
      throw std::runtime_error(
          "Target O must be finite and non-negative in build_species_plan."
      );
    }
    
    if (len == 0u) {
      P.O_n[sp] = 0.0;
    } else {
      P.O_n[sp] = static_cast<double>(len) * o;
    }

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
  
  // Fill flattened reverse pool
  P.cell_species_pool.assign(total_memberships, -1);
  std::vector<uint32_t> cursor = P.cell_species_start;
  
  for (int sp = 0; sp < n_species; ++sp) {
    const uint32_t start = P.seed_start[sp];
    const uint32_t len   = P.seed_len[sp];
    
    for (uint32_t i = 0; i < len; ++i) {
      const uint32_t cell = static_cast<uint32_t>(P.seed_pool[start + i]);
      P.cell_species_pool[cursor[cell]++] = sp;
    }
  }

  return P;
}
