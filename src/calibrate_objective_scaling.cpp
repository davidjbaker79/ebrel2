//--------------------------- Calibrate objective scaling --------------------

#include "calibrate_objective_scaling.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>

//-------------------------- Local helpers -------------------------------------

namespace {

// Helper function for getting cell index for dim_x * dim_y * n_actions vector
inline int idx(const int cell, const int action, const int n_cells) {
  return action * n_cells + cell;
}
 
 void validate_common_inputs(
     const std::vector<uint8_t>& U,
     const std::vector<int8_t>& X_reference,
     const int n_actions,
     const int dim_x,
     const int dim_y
 ) {
   if (dim_x <= 0 || dim_y <= 0) {
     throw std::runtime_error("dim_x and dim_y must both be positive.");
   }
   if (n_actions <= 0) {
     throw std::runtime_error("n_actions must be positive.");
   }
   
   const int n_cells = dim_x * dim_y;
   
   if (static_cast<int>(X_reference.size()) != n_cells) {
     throw std::runtime_error(
         "X_reference must have length dim_x * dim_y.");
   }
   if (static_cast<int>(U.size()) != n_actions * n_cells) {
     throw std::runtime_error(
         "U must have length n_actions * dim_x * dim_y.");
   }
 }

  std::vector<int> eligible_cells_for_action(
      const std::vector<uint8_t>& U,
      const std::vector<int8_t>& X_reference,
      const int action,
      const int n_cells
  ) {
    std::vector<int> eligible;
    eligible.reserve(n_cells);
    
    for (int cell = 0; cell < n_cells; ++cell) {
      // The calibration move represents an addition to the fixed reference
      // solution, rather than a replacement of an existing planned action.
      if (X_reference[cell] != static_cast<int8_t>(-1)) {
        continue;
      }
      
      if (U[idx(cell, action, n_cells)] == static_cast<uint8_t>(0)) {
        eligible.push_back(cell);
      }
    }
    
    return eligible;
  }
  
  void append_if_nonzero(
      const double value,
      const double epsilon,
      std::vector<double>& out
  ) {
    if (std::isfinite(value) && std::abs(value) > epsilon) {
      out.push_back(std::abs(value));
    }
  }
  
  // `values` is passed by value because nth_element reorders it.
  double checked_median_scale(std::vector<double> values, const double epsilon) {
    values.erase(
      std::remove_if(
        values.begin(), values.end(),
        [epsilon](const double x) {
          return !std::isfinite(x) || x <= epsilon;
        }),
        values.end());
    
    if (values.empty()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    
    const std::size_t n = values.size();
    const std::size_t middle = n / 2;
    
    std::nth_element(
      values.begin(),
      values.begin() + static_cast<std::ptrdiff_t>(middle),
      values.end());
    
    const double upper_middle = values[middle];
    
    if ((n % 2U) == 1U) {
      return upper_middle;
    }
    
    // nth_element guarantees all elements before `middle` are <= upper_middle.
    const double lower_middle = *std::max_element(
      values.begin(),
      values.begin() + static_cast<std::ptrdiff_t>(middle));
    
    return 0.5 * (lower_middle + upper_middle);
  }
  
  ActionCalibrationSummary make_action_summary(
      const int action,
      const int n_eligible_cells,
      const int n_moves_requested,
      const int n_moves_valid,
      const std::vector<double>& abs_cost,
      const std::vector<double>& abs_config,
      const std::vector<double>& abs_target,
      const double epsilon
  ) {
    ActionCalibrationSummary out;
    out.action = action;
    out.n_eligible_cells = n_eligible_cells;
    out.n_moves_requested = n_moves_requested;
    out.n_moves_valid = n_moves_valid;
    
    out.n_cost_nonzero = static_cast<int>(abs_cost.size());
    out.n_config_nonzero = static_cast<int>(abs_config.size());
    out.n_target_nonzero = static_cast<int>(abs_target.size());
    
    const double cost = checked_median_scale(abs_cost, epsilon);
    const double config = checked_median_scale(abs_config, epsilon);
    const double target = checked_median_scale(abs_target, epsilon);
    
    // NA-like diagnostics are less convenient in C++ output.  Zero means that
    // this particular action had no non-zero samples for that component.  The
    // pooled calibration below still fails if an entire component is inactive.
    out.median_abs_cost = std::isfinite(cost) ? cost : 0.0;
    out.median_abs_config = std::isfinite(config) ? config : 0.0;
    out.median_abs_target = std::isfinite(target) ? target : 0.0;
    
    return out;
  }
  
} // namespace


//-------------------------- Main functions ------------------------------------

bool sample_calibration_change(
    const std::vector<uint8_t>& U,
    const std::vector<int8_t>& X_reference,
    const int action,
    const int n_actions,
    const int dim_x,
    const int dim_y,
    std::mt19937& rng,
    CellChange& out_change
) {
  validate_common_inputs(U, X_reference, n_actions, dim_x, dim_y);
  
  if (action < 0 || action >= n_actions) {
    throw std::runtime_error("action must be in [0, n_actions).");
  }
  
  const int n_cells = dim_x * dim_y;
  const std::vector<int> eligible = eligible_cells_for_action(
    U, X_reference, action, n_cells);
  
  if (eligible.empty()) {
    return false;
  }
  
  std::uniform_int_distribution<int> draw(
      0, static_cast<int>(eligible.size()) - 1);
  
  const int cell = eligible[draw(rng)];
  out_change = CellChange{
    cell,
    X_reference[cell],
               static_cast<int8_t>(action)
  };
  
  return true;
}

ObjectiveScalingCalibration calibrate_objective_scales(
    const std::vector<uint8_t>& U,
    const std::vector<int8_t>& X_reference,
    const int n_actions,
    const int dim_x,
    const int dim_y,
    const int n_samples_per_action,
    const RawObjectiveEvaluator& evaluate_raw_components,
    std::mt19937& rng,
    const double epsilon
) {
  validate_common_inputs(U, X_reference, n_actions, dim_x, dim_y);
  
  if (n_samples_per_action <= 0) {
    throw std::runtime_error("n_samples_per_action must be positive.");
  }
  if (!evaluate_raw_components) {
    throw std::runtime_error("evaluate_raw_components must be supplied.");
  }
  if (!std::isfinite(epsilon) || epsilon < 0.0) {
    throw std::runtime_error("epsilon must be finite and non-negative.");
  }
  
  const int n_cells = dim_x * dim_y;
  
  // The base state is evaluated once.  Existing habitat and all fixed
  // species/dispersal context should already be incorporated by the evaluator.
  const RawObjectiveComponents base = evaluate_raw_components(X_reference);
  
  if (!std::isfinite(base.cost) ||
      !std::isfinite(base.config) ||
      !std::isfinite(base.target)) {
      throw std::runtime_error(
          "Raw objective evaluator returned a non-finite value for X_reference.");
  }
  
  std::vector<double> pooled_abs_cost;
  std::vector<double> pooled_abs_config;
  std::vector<double> pooled_abs_target;
  
  pooled_abs_cost.reserve(static_cast<std::size_t>(n_actions) *
    static_cast<std::size_t>(n_samples_per_action));
  pooled_abs_config.reserve(static_cast<std::size_t>(n_actions) *
    static_cast<std::size_t>(n_samples_per_action));
  pooled_abs_target.reserve(static_cast<std::size_t>(n_actions) *
    static_cast<std::size_t>(n_samples_per_action));
  
  ObjectiveScalingCalibration out;
  out.scales.n_actions_requested = n_actions;
  out.scales.n_moves_requested = n_actions * n_samples_per_action;
  out.by_action.reserve(n_actions);
  
  for (int action = 0; action < n_actions; ++action) {
    const std::vector<int> eligible = eligible_cells_for_action(
      U, X_reference, action, n_cells);
    
    std::vector<double> action_abs_cost;
    std::vector<double> action_abs_config;
    std::vector<double> action_abs_target;
    
    if (eligible.empty()) {
      out.by_action.push_back(make_action_summary(
          action, 0, n_samples_per_action, 0,
          action_abs_cost, action_abs_config, action_abs_target, epsilon));
      continue;
    }
    
    ++out.scales.n_actions_sampled;
    
    // Sample without replacement whenever possible.  It makes the geographic
    // sample broader and avoids duplicate points for small action masks.
    std::vector<int> sampled_cells = eligible;
    std::shuffle(sampled_cells.begin(), sampled_cells.end(), rng);
    
    const int n_draws = std::min(
      n_samples_per_action,
      static_cast<int>(sampled_cells.size()));
    
    for (int draw = 0; draw < n_draws; ++draw) {
      const int cell = sampled_cells[draw];
      
      std::vector<int8_t> X_proposed = X_reference;
      X_proposed[cell] = static_cast<int8_t>(action);
      
      const RawObjectiveComponents proposal =
        evaluate_raw_components(X_proposed);
      
      if (!std::isfinite(proposal.cost) ||
          !std::isfinite(proposal.config) ||
          !std::isfinite(proposal.target)) {
          std::ostringstream message;
        message << "Raw objective evaluator returned a non-finite value for "
                << "calibration action " << action << " at cell " << cell << ".";
        throw std::runtime_error(message.str());
      }
      
      const double delta_cost = proposal.cost - base.cost;
      const double delta_config = proposal.config - base.config;
      const double delta_target = proposal.target - base.target;
      
      append_if_nonzero(delta_cost, epsilon, action_abs_cost);
      append_if_nonzero(delta_config, epsilon, action_abs_config);
      append_if_nonzero(delta_target, epsilon, action_abs_target);
      
      append_if_nonzero(delta_cost, epsilon, pooled_abs_cost);
      append_if_nonzero(delta_config, epsilon, pooled_abs_config);
      append_if_nonzero(delta_target, epsilon, pooled_abs_target);
      
      ++out.scales.n_moves_valid;
      
      if (out.scales.n_moves_valid <= 10) {
        std::cout
        << "[calibration move]"
        << " action=" << action
        << " cell=" << cell
        << " d_cost=" << delta_cost
        << " d_config=" << delta_config
        << " d_target=" << delta_target
        << " base_cost=" << base.cost
        << " prop_cost=" << proposal.cost
        << " base_config=" << base.config
        << " prop_config=" << proposal.config
        << " base_target=" << base.target
        << " prop_target=" << proposal.target
        << '\n';
      }
      
    }
    
    out.by_action.push_back(make_action_summary(
        action,
        static_cast<int>(eligible.size()),
        n_samples_per_action,
        n_draws,
        action_abs_cost,
        action_abs_config,
        action_abs_target,
        epsilon));
  }
  
  out.scales.cost = checked_median_scale(pooled_abs_cost, epsilon);
  out.scales.config = checked_median_scale(pooled_abs_config, epsilon);
  out.scales.target = checked_median_scale(pooled_abs_target, epsilon);
  
  out.scales.n_cost_nonzero = static_cast<int>(pooled_abs_cost.size());
  out.scales.n_config_nonzero = static_cast<int>(pooled_abs_config.size());
  out.scales.n_target_nonzero = static_cast<int>(pooled_abs_target.size());
  
  if (!std::isfinite(out.scales.cost) ||
      !std::isfinite(out.scales.config) ||
      !std::isfinite(out.scales.target)) {
      std::ostringstream message;
    message << "Objective-scale calibration failed: no finite, non-zero "
            << "sampled change for";
    
    if (!std::isfinite(out.scales.cost)) {
      message << " cost";
    }
    if (!std::isfinite(out.scales.config)) {
      message << " configuration";
    }
    if (!std::isfinite(out.scales.target)) {
      message << " target";
    }
    message << ".";
    throw std::runtime_error(message.str());
  }
  
  out.scales.calibrated = true;
  return out;
}
