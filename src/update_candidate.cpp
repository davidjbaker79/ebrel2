//--------------------------- Update candidate ---------------------------------

#include "update_candidate.h"
#include <random>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <cmath>

//-------------------------- Local helpers -------------------------------------

namespace {

  // Helper function for getting cell index for dim_x * dim_y * n_actions vector
  static inline int idx(int cell, int a, int n_cells) {
    return a * n_cells + cell;
  }

}

//-------------------------- Main functions ------------------------------------


// Update candidate solution based on W (distance decay weighting)
CandidateUpdate update_candidate(
    const std::vector<double>& W,
    const std::vector<uint8_t>& U,
    const std::vector<int8_t>& candidate_in,
    double step_proportion,
    double step_probability,
    double improve_total_weight,
    int n_actions,
    int n_habitats,
    int dim_x,
    int dim_y,
    int rng_seed = -1
) {
  const int n_cells = dim_x * dim_y;
  
  // Check input dimensions
  if ((int)candidate_in.size() != n_cells)
    throw std::runtime_error("candidate_in must have length dim_x*dim_y");
  if ((int)W.size() != n_habitats * n_cells)
    throw std::runtime_error("W must have length n_habitats*dim_x*dim_y");
  if ((int)U.size() != n_actions * n_cells)
    throw std::runtime_error("U must have length n_actions*dim_x*dim_y");
  if (n_actions != n_habitats + 1)
    throw std::runtime_error("Expected n_actions = n_habitats + 1");
  
  // Results to CandidateUpdate
  CandidateUpdate result;
  result.X = candidate_in; 
  
  // Check on inputs
  if (step_proportion <= 0.0 || step_probability < 0.0) return result;
  if (step_probability > 1.0) throw std::runtime_error("step_probability must be <= 1");
  
  // Collect eligible cells (any habitat available, now inc. improvement)
  std::vector<int> eligible;
  eligible.reserve(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int a = 0; a < n_actions; ++a) {
      if (U[idx(i, a, n_cells)] == 0) {
        eligible.push_back(i);
        break;
      }
    }
  }
  if (eligible.empty()) return result;
  
  // How many cells to try
  int step_size = static_cast<int>(
    std::round(step_proportion * static_cast<double>(eligible.size())));
  if (step_size <= 0) return result;
  
  // Set improvement proposal weights.
  //
  // Creation weights in W are assumed to be normalised within each
  // habitat/action layer, such that each creation habitat has total
  // proposal mass approximately equal to 1 across eligible cells.
  //
  // Give improvement a total proposal mass of improve_w_default,
  // distributed equally across all cells where improvement is available.
  const int improve_action = n_actions - 1;  // zero-indexed final action
  
  int n_improvable = 0;
  for (int i = 0; i < n_cells; ++i) {
    if (U[idx(i, improve_action, n_cells)] == 0) {
      ++n_improvable;
    }
  }
  
  const double improve_w_per_cell =
    (n_improvable > 0)
    ? improve_total_weight / static_cast<double>(n_improvable)
      : 0.0;
  
  // Per-cell proposal mass used when selecting a cell to update.
  // Creation actions use their cell-specific W values; improvement
  // contributes its equal share of the total improvement proposal mass.
  std::vector<double> cell_probs(n_cells, 0.0);
  
  for (int i = 0; i < n_cells; ++i) {
    double s = 0.0;
    
    // Creation actions
    for (int a = 0; a < n_habitats; ++a) {
      if (U[idx(i, a, n_cells)] == 0) {
        s += W[idx(i, a, n_cells)];
      }
    }
    
    // Improvement action
    if (U[idx(i, improve_action, n_cells)] == 0) {
      s += improve_w_per_cell;
    }
    
    cell_probs[i] = s;
  }
  
  // Restrict to eligible & normalize
  std::vector<double> eligible_probs;
  eligible_probs.reserve(eligible.size());
  for (int i : eligible) eligible_probs.push_back(cell_probs[i]);
  
  double total = std::accumulate(eligible_probs.begin(), eligible_probs.end(), 0.0);
  if (total <= 0.0) throw std::runtime_error("update_candidate: All values in W are zero.");
  for (double& p : eligible_probs) p /= total;
  
  // RNG
  std::mt19937 gen;
  if (rng_seed >= 0) gen.seed(static_cast<uint32_t>(rng_seed));
  else gen.seed(std::random_device{}());
  
  std::vector<double> draw_probs = eligible_probs;
  std::vector<double> action_weights(n_actions, 0.0);
  std::vector<double> outcome_probs(n_actions + 1, 0.0);
  
  // Single loop — build result.X and result.changes together
  for (int s = 0; s < step_size; ++s) {
    double remaining =
      std::accumulate(draw_probs.begin(), draw_probs.end(), 0.0); // sampling without replacement
    
    if (remaining <= 0.0) break;
    
    std::discrete_distribution<> cell_dist(draw_probs.begin(), draw_probs.end());
    const int picked = cell_dist(gen);
    const int i = eligible[picked];
    
    draw_probs[picked] = 0.0;
    
    double sum_w = 0.0;
    
    // Creation actions
    for (int a = 0; a < n_habitats; ++a) {
      const double w = (U[idx(i, a, n_cells)] == 0) ? W[idx(i, a, n_cells)] : 0.0;
      action_weights[a] = w;
      sum_w += w;
    }
    
    // Improvement action
    action_weights[improve_action] =
      (U[idx(i, improve_action, n_cells)] == 0) ? improve_w_per_cell : 0.0;
    sum_w += action_weights[improve_action];
    
    outcome_probs[0] = 1.0 - step_probability;
    
    if (sum_w > 0.0) {
      const double scale = step_probability / sum_w;
      for (int a = 0; a < n_actions; ++a) {
        outcome_probs[a + 1] = action_weights[a] * scale;
      }
    } else {
      std::fill(outcome_probs.begin() + 1, outcome_probs.end(), 0.0);
      outcome_probs[0] = 1.0;
    }
    
    std::discrete_distribution<> outcome_dist(outcome_probs.begin(), outcome_probs.end());
    const int outcome = outcome_dist(gen);
    
    const int8_t new_action =
      (outcome == 0) ? static_cast<int8_t>(-1)
        : static_cast<int8_t>(outcome - 1);
    
    if (new_action != result.X[i]) {
      result.changes.push_back({i, result.X[i], new_action});
      result.X[i] = new_action;
    }
  }
  
  return result;
}