//-------------------------- Run Ebrel model -----------------------------------

#include "run_ebrel.h"
#include "simulated_annealing.h"
#include "objective_utils.h"
#include "species_plan.h"
#include "dispersal_utils.h"
#include "generate_x_zero.h"
#include "optimisation_utils.h"


#include <numeric>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdint>

//---------------------- Local helper functions --------------------------------

namespace {

inline void ensure(bool ok, const char* msg) {
  if (!ok) {
    throw std::invalid_argument(msg);
  }
}
  
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) *
      static_cast<std::size_t>(dim_y);
  }
  
  inline void validate_objective_settings(const RunEBRELInput& in) {
    
    // User-defined relative preferences.
    ensure(std::isfinite(in.alpha), "alpha must be finite");
    ensure(std::isfinite(in.beta),  "beta must be finite");
    ensure(std::isfinite(in.gamma), "gamma must be finite");
    
    ensure(in.alpha >= 0.0, "alpha must be non-negative");
    ensure(in.beta  >= 0.0, "beta must be non-negative");
    ensure(in.gamma >= 0.0, "gamma must be non-negative");
    
    ensure(
      (in.alpha > 0.0) || (in.beta > 0.0) || (in.gamma > 0.0),
      "At least one of alpha, beta, or gamma must be positive"
    );
    
    // Internally calibrated raw-component scales.
    ensure(
      std::isfinite(in.objective_scales.cost) &&
        in.objective_scales.cost > 0.0,
        "objective_scales.cost must be finite and positive"
    );
    
    ensure(
      std::isfinite(in.objective_scales.config) &&
        in.objective_scales.config > 0.0,
        "objective_scales.config must be finite and positive"
    );
    
    ensure(
      std::isfinite(in.objective_scales.target) &&
        in.objective_scales.target > 0.0,
        "objective_scales.target must be finite and positive"
    );
  }
  
  inline void validate_shapes(const RunEBRELInput& in) {
    
    const std::size_t cells = n_cells(in.dim_x, in.dim_y);
    
    const std::size_t Asz =
      cells * static_cast<std::size_t>(in.n_actions);
    
    const std::size_t Ssz =
      cells * static_cast<std::size_t>(in.n_species);
    
    const std::size_t SxHsz =
      static_cast<std::size_t>(in.n_habitats) *
      static_cast<std::size_t>(in.n_species);
    
    ensure(
      in.dim_x > 0 && in.dim_y > 0,
      "dim_x and dim_y must be positive"
    );
    
    ensure(in.n_actions  > 0, "n_actions must be positive");
    ensure(in.n_species  > 0, "n_species must be positive");
    ensure(in.n_habitats > 0, "n_habitats must be positive");
    
    ensure(
      in.U.size() == Asz,
      "U size mismatch (expected dim_x*dim_y*n_actions)"
    );
    
    ensure(
      in.C.size() == Asz,
      "C size mismatch (expected dim_x*dim_y*n_actions)"
    );
    
    ensure(
      in.E.size() == cells,
      "E size mismatch (expected dim_x*dim_y)"
    );
    
    ensure(
      in.O.size() == static_cast<std::size_t>(in.n_species),
      "O size mismatch (expected n_species)"
    );
    
    ensure(
      in.SD.size() == Ssz,
      "SD size mismatch (expected dim_x*dim_y*n_species)"
    );
    
    ensure(
      in.SxH.size() == SxHsz,
      "SxH size mismatch (expected n_habitats*n_species)"
    );
    
    ensure(
      in.D.size() == static_cast<std::size_t>(in.n_species),
      "D size mismatch (expected n_species)"
    );
    
    if (!in.X0.empty()) {
      ensure(
        in.X0.size() == cells,
        "X0 size mismatch (expected dim_x*dim_y)"
      );
    }
    
    validate_objective_settings(in);
  }
  
} // namespace


//--------------------- Main functions -----------------------------------------

RunEBRELResult run_ebrel(const RunEBRELInput& in, 
                         const RunEBRELOptions& opt) {
  validate_shapes(in);

  if (opt.verbose) {
    std::cout << "[run_ebrel] Starting SA with "
              << in.dim_x << "x" << in.dim_y
              << ", n_actions =" << in.n_actions 
              << ", n_species =" << in.n_species
              << ", iterations=" << opt.n_iterations 
              << std::endl;
  }
  
  SAResult sa = simulated_annealing(
    in.X0, in.W, in.U, in.C, in.E,
    in.Etiles_per_h,
    in.cell_r, in.cell_c,
    in.rowruns_cache,
    in.species_plan,
    in.improve_count,
    in.n_actions, in.n_species, in.n_habitats,
    in.dim_x, in.dim_y,
    in.universal_disp_thres,
    in.max_disp_steps,
    in.roi_cap,
    in.LM,
    in.row_first_land, in.row_last_land, in.col_first_land, in.col_last_land,
    in.alpha, in.beta, in.gamma,
    in.objective_scales.cost, in.objective_scales.config, in.objective_scales.target,
    in.improve_w_default,
    opt.step_proportion,
    opt.step_probability,
    opt.n_iterations,
    opt.temp,
    opt.cooling_rate_c,
    opt.lam_enabled,
    opt.lam_target_mid,
    opt.lam_target_final,
    opt.lam_hold_frac,
    opt.lam_p,
    opt.min_iterations,
    opt.acceptance_window,
    opt.acceptance_thres,
    opt.iter_no_improve,
    opt.improve_eps,
    opt.write_every,
    opt.trace_file,
    opt.verbose
  );

  RunEBRELResult out;
  out.X_best         = std::move(sa.X_best);
  out.H_best         = sa.H_best;
  out.H_trace        = std::move(sa.H_trace);
  out.F_trace        = std::move(sa.F_trace);
  out.F1_trace       = std::move(sa.F1_trace);
  out.F2_trace       = std::move(sa.F2_trace);
  out.iterations_run = static_cast<int>(out.H_trace.size());
  out.g_best         = std::move(sa.g_best);
  out.g_create_best  = std::move(sa.g_create_best);
  out.g_improve_best  = std::move(sa.g_improve_best);

  // --- NEW: diagnostics ---
  out.acc_rate_trace   = std::move(sa.diag.acceptance_rate_trace);
  out.early_stop_iter  = sa.diag.early_stop_iter;

  // Prefer SA’s attempted_total; fall back to trace length if zero
  out.proposals = (sa.diag.attempted_total > 0)
    ? sa.diag.attempted_total
  : static_cast<int>(out.H_trace.size());
  out.accepted = sa.diag.accepted_total;
  out.overall_acc = (out.proposals > 0)
    ? static_cast<double>(out.accepted) / out.proposals
  : std::numeric_limits<double>::quiet_NaN();

  // Timings
  out.iter_ms_total = sa.diag.iter_ms_total;
  out.iter_count    = sa.diag.iter_count;

  if (opt.verbose) {
    std::cout << "[run_ebrel] Done. Best H = " << out.H_best
              << " (cand evals recorded: " << out.iterations_run << ")"
              << std::endl;

  }
  return out;
}
