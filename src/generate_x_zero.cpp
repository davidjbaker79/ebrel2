//------------------------- X0 functions ---------------------------------------

#include "generate_x_zero.h"

#include <vector>
#include <random>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cctype>
#include <cstdint>

//-------------------------- Local helper functions ----------------------------

namespace {

  // Number of cells
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  // For safety check in [0,1]
  inline bool in_unit_interval(double p) {
    return std::isfinite(p) && p >= 0.0 && p <= 1.0;
  }

}

//----------------------------- Main functions ---------------------------------

/* Naive initialize of habitat configuration matrix X0 - includes options for 
   habitat creation and improvement */
std::vector<int8_t> generate_X0_CI(const std::vector<uint8_t>& U,
                                   int n_actions,
                                   int dim_x,
                                   int dim_y,
                                   double base_prob,
                                   int rng_seed) {
  // --- validation ---
  if (dim_x <= 0 || dim_y <= 0 || n_actions <= 0) {
    throw std::invalid_argument("dim_x, dim_y, and n_actions must be positive.");
  }
  
  if (!in_unit_interval(base_prob)) {
    throw std::invalid_argument("base_prob must be in [0, 1].");
  }
  
  const std::size_t cells    = n_cells(dim_x, dim_y);
  const std::size_t expected = cells * static_cast<std::size_t>(n_actions);
  if (U.size() != expected) {
    throw std::invalid_argument("U size mismatch: got " + std::to_string(U.size()) +
                                ", expected " + std::to_string(expected));
  }

  // If base_prob == 1 set everything to everything empty (-1)
  std::vector<int8_t> X0(cells, static_cast<int8_t>(-1));
  if (base_prob >= 1.0) return X0;

  // RNG: deterministic if rng_seed >= 0, random_device otherwise
  std::mt19937 rng;
  if (rng_seed >= 0) {
    rng.seed(static_cast<std::mt19937::result_type>(rng_seed));
  } else { 
    std::random_device rd; rng.seed(rd()); 
  }

  std::uniform_real_distribution<double> unif01(0.0, 1.0);
  std::uniform_int_distribution<int> habitat_dist(0, (n_actions-1)); // This is the n_actions + improvement (zero indexed)

  // For each tile: leave empty with probability base_prob; else assign an available habitat
  for (std::size_t tile = 0; tile < cells; ++tile) {
    if (unif01(rng) < base_prob) continue;   // empty

    const int h0 = habitat_dist(rng);

    // Check habitats in rotated order, pick first available
    for (int k = 0; k < (n_actions-1); ++k) {
      const int h = (h0 + k) % (n_actions-1);
      const std::size_t idx = static_cast<std::size_t>(h) * cells + tile; // U is habitat-major
      if (U[idx] == 0u) {
        X0[tile] = static_cast<int8_t>(h);
        break;
      }
    }
    // if none available, remains -1
  }

  return X0;
}
