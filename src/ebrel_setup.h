//--------------------------- Ebrel setup functions ----------------------------

// These functions are for creating an Ebrel object with validation. The Ebrel
// object is then supplied to the run_ebrel function.

// Memory layout note for reference:
// - Action-major 3D fields (E, C, P): idx = h * (dim_x*dim_y) + cell
// - Species-major 3D field (SD):     idx = s * (dim_x*dim_y) + cell
// - SxH (n_h * n_s):                 idx = h * n_s + s
// - O is species-level:              size == n_s

#ifndef EBREL_SETUP_H
#define EBREL_SETUP_H

#include <vector>
#include <stdexcept>
#include <string>
#include <cstddef>
#include <cmath>

#include "offsets.h"

// Validate input dimensions for creating the Ebrel class object
inline void validate_dimensions(const std::vector<int8_t>& E,
                                const std::vector<double>& C,
                                const std::vector<double>& SD,
                                const std::vector<int>& D,
                                const std::vector<double>& SxH,
                                const std::vector<double>& O,
                                const std::vector<float>& action_weight,
                                const std::vector<uint8_t>& LM,
                                int dim_x,
                                int dim_y,
                                int n_actions,
                                int n_species) {
  if (dim_x <= 0 || dim_y <= 0 || n_actions <= 0 || n_species <= 0) {
    throw std::runtime_error("dim_x, dim_y, n_actions, and n_species must all be positive.");
  }
  
  // Define number of habitat creation options (last action is improvement)
  int n_habitats = n_actions - 1;

  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const std::size_t Asz   = cells * static_cast<std::size_t>(n_actions);
  const std::size_t Ssz   = cells * static_cast<std::size_t>(n_species);
  const std::size_t SxHsz = static_cast<std::size_t>(n_habitats) * static_cast<std::size_t>(n_species);

  if (E.size()  != cells)
    throw std::runtime_error("E size mismatch: got " + std::to_string(E.size())  +
                             ", expected " + std::to_string(cells));
  if (C.size()  != Asz)
    throw std::runtime_error("C size mismatch: got " + std::to_string(C.size())  +
                             ", expected " + std::to_string(Asz) + " (one per cell).");
  if (SD.size() != Ssz)
    throw std::runtime_error("SD size mismatch: got " + std::to_string(SD.size()) +
                             ", expected " + std::to_string(Ssz));
  if (D.size()  != static_cast<std::size_t>(n_species))
    throw std::runtime_error("D size mismatch: got " + std::to_string(D.size())  +
                             ", expected " + std::to_string(n_species));
  if (SxH.size()!= SxHsz)
    throw std::runtime_error("SxH size mismatch: got " + std::to_string(SxH.size()) +
                             ", expected " + std::to_string(SxHsz));
  if (O.size()  != static_cast<std::size_t>(n_species))
    throw std::runtime_error("O size mismatch: got " + std::to_string(O.size()) +
                             ", expected " + std::to_string(n_species) + " (species-level).");
  if (LM.size() != cells)
    throw std::runtime_error("LM size mismatch: got " + std::to_string(LM.size()) +
                             ", expected " + std::to_string(cells) + " (one per cell).");
  if (action_weight.size() != static_cast<std::size_t>(n_actions))
    throw std::runtime_error("Action_weight must have length n_actions");
  
  // In validate_dimensions — explicit E range check
  for (std::size_t i = 0; i < E.size(); ++i) {
    if (E[i] != -1 && E[i] >= n_habitats) {
      throw std::runtime_error(
          "E contains value " + std::to_string(E[i]) +
            " at cell " + std::to_string(i) +
            " which is >= n_habitats (" + std::to_string(n_habitats) + ")");
    }
  }
  
}

// Create unavailable mask from C
[[nodiscard]] std::vector<uint8_t> get_unavailable_mask(const std::vector<double>& C,
                                                       int n_actions,
                                                       int dim_x,
                                                       int dim_y,
                                                       double sentinel = 1e10);

// Compute land spans by row and col for skipping sea
inline void compute_land_spans(const std::vector<uint8_t>& LM,
                               int dim_x,
                               int dim_y,
                               std::vector<int16_t>& row_first_land,
                               std::vector<int16_t>& row_last_land,
                               std::vector<int16_t>& col_first_land,
                               std::vector<int16_t>& col_last_land)
{
  // One entry per row / column
  row_first_land.assign(dim_y, dim_x);  // "no land" => first = dim_x
  row_last_land.assign(dim_y, -1);      // "no land" => last = -1

  col_first_land.assign(dim_x, dim_y);  // "no land" => first = dim_y
  col_last_land.assign(dim_x, -1);      // "no land" => last = -1

  // Scan LM row-major: g = r * dim_x + c
  for (int r = 0; r < dim_y; ++r) {
    std::size_t base = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x);
    for (int c = 0; c < dim_x; ++c) {
      std::size_t g = base + static_cast<std::size_t>(c);
      if (!LM[g]) continue;  // sea

      if (c < row_first_land[r]) row_first_land[r] = c;
      if (c > row_last_land[r])  row_last_land[r]  = c;

      if (r < col_first_land[c]) col_first_land[c] = r;
      if (r > col_last_land[c])  col_last_land[c]  = r;
    }
  }
}

// Build Ebrel class object outputs
inline void create_ebrel_class_object(const std::vector<int8_t>& E,
                                      const std::vector<double>& C,
                                      const std::vector<double>& SD,
                                      const std::vector<int>&    D,
                                      const std::vector<double>& SxH,
                                      const std::vector<double>& O,
                                      const std::vector<float>& action_weight, 
                                      const std::vector<uint8_t>& LM,
                                      int universal_disp_thres,
                                      int dim_x,
                                      int dim_y,
                                      int n_actions,
                                      int n_species,
                                      int n_habitats,
                                      double sentinel,
                                      double sigma_in,
                                      double& sigma_out,
                                      std::vector<uint8_t>& U_out,
                                      std::vector<int16_t>& row_first_land,
                                      std::vector<int16_t>& row_last_land,
                                      std::vector<int16_t>& col_first_land,
                                      std::vector<int16_t>& col_last_land,
                                      std::vector<int>& cell_r_out,
                                      std::vector<int>& cell_c_out,
                                      std::vector<std::vector<std::size_t>>& Etiles_per_h_out,
                                      RowRunsCache& rowruns_cache_out)
{
  // ---- Check sizes, including LM ----
  validate_dimensions(E, C, SD, D, SxH, O, action_weight, LM, dim_x, dim_y, n_actions, n_species);

  // ---- Derive land spans from the LM mask ----
  compute_land_spans(LM, dim_x, dim_y,
                     row_first_land, row_last_land,
                     col_first_land, col_last_land);

  // ---- Existing unavailable-mask logic ----
  U_out = get_unavailable_mask(C, n_actions, dim_x, dim_y, sentinel);

  // ---- Apply land mask to U ----
  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  for (int a = 0; a < n_actions; ++a) {
    const std::size_t off = static_cast<std::size_t>(a) * cells;
    for (std::size_t cell = 0; cell < cells; ++cell) {
      if (LM[cell] == 0) U_out[off + cell] = 1;
    }
  }

  // ---- Check that no cells outside LM have E ----
  for (std::size_t cell = 0; cell < cells; ++cell) {
    if (LM[cell] == 0 && E[cell] != -1) {
      throw std::runtime_error("E must be -1 where LM==0 (cell=" + std::to_string(cell) + ")");
    }
  }

  // ---- Calculate sigma = 1 / mean(D>0) ----
  double sigma;

  if (sigma_in > 0.0 && std::isfinite(sigma_in)) {
    sigma = sigma_in;
  } else {
    if (D.empty()) {
      throw std::runtime_error("D must not be empty when computing default sigma");
    }
    double sumD = 0.0;
    int cnt = 0;
    for (int d : D) {
      if (d > 0) { sumD += static_cast<double>(d); ++cnt; }
    }
    if (cnt == 0) {
      throw std::runtime_error("All values in D are <= 0; cannot compute default sigma");
    }

    const double meanD = sumD / static_cast<double>(cnt);
    sigma = 1.0 / meanD;

    if (!(sigma > 0.0) || !std::isfinite(sigma)) {
      throw std::runtime_error("Computed default sigma is not positive/finite");
    }
  }
  sigma_out = sigma;

  // ---- Precompute row/col coords [TESTED in codePad & WORKS] ----
  cell_r_out.resize(cells);
  cell_c_out.resize(cells);
  for (std::size_t k = 0; k < cells; ++k) {
    cell_r_out[k] = static_cast<int>(k / static_cast<std::size_t>(dim_x));
    cell_c_out[k] = static_cast<int>(k % static_cast<std::size_t>(dim_x));
  }

  // ---- Pre-compute Etiles for compute_F2 ----
  Etiles_per_h_out.assign(static_cast<std::size_t>(n_habitats), {});
  for (std::size_t tile = 0; tile < cells; ++tile) {
    const int h = E[tile];
    if (h >= 0) {
      if (h >= n_habitats) {
        throw std::runtime_error("E contains habitat index >= n_actions - 1");
      }
      Etiles_per_h_out[static_cast<std::size_t>(h)].push_back(tile);
    }
  }

  // ---- Build offsets cache ----
  OffsetsCache offsets_cache = build_offsets_cache_round(universal_disp_thres);

  // ---- Build row runs cache from offsets cache ----
  rowruns_cache_out = build_rowruns_cache_from_offsets(offsets_cache);

}

#endif  // EBREL_SETUP_H
