//--------------------------- Ebrel setup functions ----------------------------

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstddef>
#include <cstdint>


// ------------------------- Local helper functions ----------------------------

namespace {

  // Just to ensure that unavailable areas aren't misclassified due to small numerical difference
  inline bool approx_equal(double x, double y, double rel = 1e-12, double abs_tol = 1e-12) {
    const double diff  = std::abs(x - y);
    const double scale = std::max(std::abs(x), std::abs(y));
    return diff <= std::max(abs_tol, rel * scale);
  }

}

//---------------------- Main functions ----------------------------------------

// Create binary mask: 1 if cost ≈ sentinel, else 0
std::vector<uint8_t> get_unavailable_mask(const std::vector<double>& C,
                                          int n_actions,
                                          int dim_x,
                                          int dim_y,
                                          double sentinel) {
  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const std::size_t expected = cells * static_cast<std::size_t>(n_actions);
  if (C.size() != expected) {
    throw std::runtime_error("C size mismatch in get_unavailable_mask: got " +
                             std::to_string(C.size()) + ", expected " + std::to_string(expected));
  }

  std::vector<uint8_t> U(expected);
  for (std::size_t i = 0; i < expected; ++i) {
    U[i] = approx_equal(C[i], sentinel) ? uint8_t{1} : uint8_t{0};
  }
  return U;
}

