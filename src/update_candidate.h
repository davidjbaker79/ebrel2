//--------------------------- Update candidate ---------------------------------

#ifndef UPDATE_CANDIDATE_H
#define UPDATE_CANDIDATE_H

#include <vector>
#include <cstdint>

/**
 * Update candidate solution for simulated annealing
 *
 * All matrices are flattened arrays:
 * - W: per-(cell, habitat) weights, typically column-normalized.
 * - U: availability mask (1 = unavailable, 0 = available).
 * - candidate_in: length n_cells. Each element is -1 for no action, or an action index in [0, n_actions).
 *
 * Returns a new flattened candidate vector after applying updates.
 */

// Return type captures both the new X and the cells that changed
struct CellChange {
  int     cell;
  int8_t  old_action;
  int8_t  new_action;
};

struct CandidateUpdate {
  std::vector<int8_t>    X;
  std::vector<CellChange> changes;
};

// Update candidate
CandidateUpdate update_candidate(
    const std::vector<double>& W, // [n_cells * n_habitats]
    const std::vector<uint8_t>& U, // [n_cells * n_actions], 1 = unavailable
    const std::vector<int8_t>& candidate_in,  // length n_cells, one-hot
    double step_proportion, // proportion of eligible cells to update
    double step_probability, // probability of assigning any actions
    double improve_total_weight, // default value for the total weight of selection prob for improvement     
    int n_actions,
    int n_habitats,
    int dim_x,
    int dim_y,
    int rng_seed);

#endif // UPDATE_CANDIDATE_H
