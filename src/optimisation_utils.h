//--------------------- Compute distance weights -------------------------------

#ifndef OPTIMISATION_UTILS_H
#define OPTIMISATION_UTILS_H

#include <vector>
#include <chrono>
#include <stdexcept>

#include "species_plan.h"

/**
 * Compute distance–decay weights W for each (cell, habitat).
 *
 * Layout:
 * - Grid has dim_y rows and dim_x cols; n_cells = dim_x * dim_y.
 * - Input U is  flattened in **habitat-major** order:
 *     index(h, cell) = h * n_cells + cell
 *   i.e., n_habitats contiguous blocks of length n_cells (one block per habitat).
 * - Output W uses the same habitat-major layout and size (n_habitats * n_cells).
 *
 * Distances:
 * - Multi-source BFS from all seed cells where E==1 for a given habitat.
 * - 8-connected neighbourhood with **unit-cost** steps (using Chebyshev).
 *
 * Weights:
 * - For each cell: w = exp(-sigma * d), where d is the integer hop distance.
 * - Cells are masked to zero weight if U==1 (unavailable) or E==1 (already present).
 * - For each habitat h, weights W[h,*] are normalized to sum to 1 if any positive mass exists.
 *
 */
std::vector<double> compute_distance_weights(
    const std::vector<int8_t>& E,  // size: n_cells, values in [-1, 0..n_habitats-1]
    const std::vector<uint8_t>& U,  // size: n_actions * n_cells (habitat-major) but ignore the improvements
    int n_habitats,
    int dim_x,
    int dim_y,
    double sigma
);

// For benchmarking
inline double ms_since(const std::chrono::steady_clock::time_point& t0) {
  using namespace std::chrono;
  return duration_cast<duration<double, std::milli>>(steady_clock::now() - t0).count();
}


// Convert dense action-major X (action*cells + g) to per-cell action id
// Rule: argmax per cell; if max <= eps → -1
inline std::vector<int> dense_X_to_cell_action(const std::vector<double>& X,
                                                int n_actions,
                                                std::size_t cells,
                                                double eps = 0.0)
{
  const std::size_t expected =
    static_cast<std::size_t>(n_actions) * cells;

  if (X.size() != expected) {
    throw std::runtime_error(
        "`X` must have length n_actions * (dim_x * dim_y).");
  }

  std::vector<int> X_h_of_cell(cells, -1);
  const double* Xptr = X.data();

  for (std::size_t g = 0; g < cells; ++g) {
    int best_a = -1;
    double best_v = eps;

    for (int a = 0; a < n_actions; ++a) {
      const double v = Xptr[static_cast<std::size_t>(a) * cells + g];
      if (v > best_v) {
        best_v = v;
        best_a = a;
      }
    }

    X_h_of_cell[g] = best_a; // -1 if none
  }

  return X_h_of_cell;
}

// Initialise improve count helper
[[nodiscard]] std::vector<int> initialise_improve_count(
    const std::vector<int8_t>& X,
    const SpeciesPlan& P
);


#endif // OPTIMISATION_UTILS_H
