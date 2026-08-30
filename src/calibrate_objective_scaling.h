//---------------------- Calibrate objective scaling --------------------------
//
// Draws broadly distributed, one-cell additions for each action and uses the
// raw objective-component changes to derive fixed scales for cost,
// configuration, and biodiversity-target performance.
//
// The scales are internal unit conversions. They are intentionally separate
// from user preference weights (alpha, beta, gamma).

#ifndef CALIBRATE_OBJECTIVE_SCALING_H
#define CALIBRATE_OBJECTIVE_SCALING_H

#include "update_candidate.h"  // CellChange

#include <cstdint>
#include <functional>
#include <random>
#include <vector>

// Raw, unweighted and unscaled objective components.  The calibration code
// never assumes whether a component is a benefit or a penalty: it uses the
// absolute change between the reference and proposed states.
struct RawObjectiveComponents {
  double cost   = 0.0;
  double config = 0.0;
  double target = 0.0;
};

// Fixed component scales estimated during pre-SA calibration.
struct ObjectiveScales {
  double cost   = 1.0;
  double config = 1.0;
  double target = 1.0;
  
  int n_actions_requested = 0;
  int n_actions_sampled   = 0;
  int n_moves_requested   = 0;
  int n_moves_valid       = 0;
  
  int n_cost_nonzero      = 0;
  int n_config_nonzero    = 0;
  int n_target_nonzero    = 0;
  
  bool calibrated         = false;
};

// Optional diagnostics retained separately for each action.  These are not
// used in scaling; the final scales are pooled across all actions.
struct ActionCalibrationSummary {
  int action              = -1;
  int n_eligible_cells    = 0;
  int n_moves_requested   = 0;
  int n_moves_valid       = 0;
  
  int n_cost_nonzero      = 0;
  int n_config_nonzero    = 0;
  int n_target_nonzero    = 0;
  
  double median_abs_cost   = 0.0;
  double median_abs_config = 0.0;
  double median_abs_target = 0.0;
};

struct ObjectiveScalingCalibration {
  ObjectiveScales scales;
  std::vector<ActionCalibrationSummary> by_action;
};

// Evaluates raw components for a complete candidate state.  Existing habitat,
// species distributions, target accounting, and any other fixed landscape
// context should be captured by the caller's lambda / functor.
using RawObjectiveEvaluator = std::function<RawObjectiveComponents(
  const std::vector<int8_t>& X
)>;

/**
 * Randomly select one feasible, currently unselected cell for `action`.
 *
 * Availability U is flattened action-major: U[action * n_cells + cell], with
 * 0 = available and 1 = unavailable.  X_reference has length n_cells, where
 * -1 denotes an unselected planning cell and non-negative values denote an
 * assigned action.
 *
 * Returns false when there is no feasible unselected cell for the action.
 */
bool sample_calibration_change(
    const std::vector<uint8_t>& U,
    const std::vector<int8_t>& X_reference,
    int action,
    int n_actions,
    int dim_x,
    int dim_y,
    std::mt19937& rng,
    CellChange& out_change
);

/**
 * Estimate objective-component scales from one-cell, action-stratified
 * additions to a fixed reference solution.
 *
 * For every action, the function samples up to n_samples_per_action feasible
 * cells uniformly from the landscape.  Each sample is evaluated as a one-cell
 * addition to X_reference.  The returned scale for each component is the
 * pooled median of non-zero absolute component changes across all actions.
 *
 * The evaluator must return raw, unweighted, unscaled component values.
 * Throws if any component has no finite, non-zero change in the calibration
 * sample, because dividing by an arbitrary fallback would hide a modelling or
 * implementation problem.
 */
ObjectiveScalingCalibration calibrate_objective_scales(
    const std::vector<uint8_t>& U,
    const std::vector<int8_t>& X_reference,
    int n_actions,
    int dim_x,
    int dim_y,
    int n_samples_per_action,
    const RawObjectiveEvaluator& evaluate_raw_components,
    std::mt19937& rng,
    double epsilon = 1e-12
);

#endif // CALIBRATE_OBJECTIVE_SCALING_H
