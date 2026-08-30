//----------------------- R wrappers for cpp functions -------------------------

#include <Rcpp.h>
#include "ebrel_setup.h"
#include "generate_x_zero.h"
#include "objective_utils.h"
#include "dispersal_utils.h"
#include "optimisation_utils.h"
#include "update_candidate.h"
#include "simulated_annealing.h"
#include "sa_utils.h"
#include "run_ebrel.h"
#include "ebrel_builder.h"

#include <algorithm>
#include <cctype>
#include <numeric>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdint>
#include <queue>
#include <string>
#include <cstring>
#include <iostream>

using namespace Rcpp;

//-------------------------- Local helpers -------------------------------------

namespace {

inline int opt_int(const Rcpp::List& opt, const char* name, const int def) {
  return opt.containsElementNamed(name) ? Rcpp::as<int>(opt[name]) : def;
}

inline double opt_double(
    const Rcpp::List& opt,
    const char* name,
    const double def
) {
  if (!opt.containsElementNamed(name)) return def;
  
  const Rcpp::NumericVector value = opt[name];
  if (value.size() == 0 || Rcpp::NumericVector::is_na(value[0])) return def;
  
  return static_cast<double>(value[0]);
}

inline bool opt_bool(const Rcpp::List& opt, const char* name, const bool def) {
  return opt.containsElementNamed(name) ? Rcpp::as<bool>(opt[name]) : def;
}

inline int action_cell_index(const int cell, const int action, const int n_cells) {
  return action * n_cells + cell;
}

void apply_run_options(const Rcpp::List& opt, RunEBRELOptions& options) {
  if (opt.containsElementNamed("base_prob_X0"))
    options.base_prob_X0 = Rcpp::as<double>(opt["base_prob_X0"]);
  if (opt.containsElementNamed("step_proportion"))
    options.step_proportion = Rcpp::as<double>(opt["step_proportion"]);
  if (opt.containsElementNamed("step_probability"))
    options.step_probability = Rcpp::as<double>(opt["step_probability"]);
  if (opt.containsElementNamed("n_iterations"))
    options.n_iterations = Rcpp::as<int>(opt["n_iterations"]);
  if (opt.containsElementNamed("temp"))
    options.temp = Rcpp::as<double>(opt["temp"]);
  if (opt.containsElementNamed("cooling_rate_c"))
    options.cooling_rate_c = Rcpp::as<double>(opt["cooling_rate_c"]);
  if (opt.containsElementNamed("lam_enabled"))
    options.lam_enabled = Rcpp::as<bool>(opt["lam_enabled"]);
  if (opt.containsElementNamed("lam_target_mid"))
    options.lam_target_mid = Rcpp::as<double>(opt["lam_target_mid"]);
  if (opt.containsElementNamed("lam_target_final"))
    options.lam_target_final = Rcpp::as<double>(opt["lam_target_final"]);
  if (opt.containsElementNamed("lam_hold_frac"))
    options.lam_hold_frac = Rcpp::as<double>(opt["lam_hold_frac"]);
  if (opt.containsElementNamed("lam_p"))
    options.lam_p = Rcpp::as<double>(opt["lam_p"]);
  if (opt.containsElementNamed("min_iterations"))
    options.min_iterations = Rcpp::as<int>(opt["min_iterations"]);
  if (opt.containsElementNamed("acceptance_window"))
    options.acceptance_window = Rcpp::as<int>(opt["acceptance_window"]);
  if (opt.containsElementNamed("acceptance_thres"))
    options.acceptance_thres = Rcpp::as<double>(opt["acceptance_thres"]);
  if (opt.containsElementNamed("iter_no_improve"))
    options.iter_no_improve = Rcpp::as<int>(opt["iter_no_improve"]);
  if (opt.containsElementNamed("improve_eps"))
    options.improve_eps = Rcpp::as<double>(opt["improve_eps"]);
  if (opt.containsElementNamed("rng_seed"))
    options.rng_seed = Rcpp::as<int>(opt["rng_seed"]);
  if (opt.containsElementNamed("write_every"))
    options.write_every = Rcpp::as<int>(opt["write_every"]);
  if (opt.containsElementNamed("trace_file"))
    options.trace_file = Rcpp::as<std::string>(opt["trace_file"]);
  if (opt.containsElementNamed("verbose"))
    options.verbose = Rcpp::as<bool>(opt["verbose"]);
}

struct TempTuneOptions {
  bool enabled = false;
  int num_samples = 400;
  double chi0 = 0.8;
  double p = 2.0;
  double tol_logchi = 1e-3;
  int max_iters = 50;
  double T1 = -1.0; // <= 0 means heuristic initial bracket.
  int max_tries_factor = 50;
  int rng_seed = -1;
  bool verbose = false;
};

TempTuneOptions parse_temp_tune_opts(const Rcpp::List& opt) {
  TempTuneOptions out;
  out.enabled = opt_bool(opt, "tune_temp", false);
  if (!out.enabled) return out;
  
  out.num_samples = opt_int(opt, "num_samples", out.num_samples);
  out.chi0 = opt_double(opt, "chi0", out.chi0);
  out.p = opt_double(opt, "p", out.p);
  out.tol_logchi = opt_double(opt, "tol_logchi", out.tol_logchi);
  out.max_iters = opt_int(opt, "max_iters", out.max_iters);
  out.max_tries_factor = opt_int(
    opt, "max_tries_factor", out.max_tries_factor);
  out.T1 = opt_double(opt, "T1", out.T1);
  
  if (opt.containsElementNamed("seed")) {
    out.rng_seed = Rcpp::as<int>(opt["seed"]);
  } else if (opt.containsElementNamed("rng_seed")) {
    out.rng_seed = Rcpp::as<int>(opt["rng_seed"]);
  }
  
  out.verbose = opt_bool(opt, "verbose", false);
  return out;
}

std::mt19937 make_rng(const int seed) {
  if (seed >= 0) {
    return std::mt19937(static_cast<std::mt19937::result_type>(seed));
  }
  
  std::random_device rd;
  return std::mt19937(rd());
}

Rcpp::List as_r_calibration(const ObjectiveScalingCalibration& calibration) {
  Rcpp::IntegerVector action(calibration.by_action.size());
  Rcpp::IntegerVector n_eligible(calibration.by_action.size());
  Rcpp::IntegerVector n_requested(calibration.by_action.size());
  Rcpp::IntegerVector n_valid(calibration.by_action.size());
  Rcpp::IntegerVector n_cost_nonzero(calibration.by_action.size());
  Rcpp::IntegerVector n_config_nonzero(calibration.by_action.size());
  Rcpp::IntegerVector n_target_nonzero(calibration.by_action.size());
  Rcpp::NumericVector median_abs_cost(calibration.by_action.size());
  Rcpp::NumericVector median_abs_config(calibration.by_action.size());
  Rcpp::NumericVector median_abs_target(calibration.by_action.size());
  
  for (std::size_t i = 0; i < calibration.by_action.size(); ++i) {
    const ActionCalibrationSummary& x = calibration.by_action[i];
    action[i] = x.action;
    n_eligible[i] = x.n_eligible_cells;
    n_requested[i] = x.n_moves_requested;
    n_valid[i] = x.n_moves_valid;
    n_cost_nonzero[i] = x.n_cost_nonzero;
    n_config_nonzero[i] = x.n_config_nonzero;
    n_target_nonzero[i] = x.n_target_nonzero;
    median_abs_cost[i] = x.median_abs_cost;
    median_abs_config[i] = x.median_abs_config;
    median_abs_target[i] = x.median_abs_target;
  }
  
  const ObjectiveScales& s = calibration.scales;
  
  return Rcpp::List::create(
    Rcpp::Named("scale_cost") = s.cost,
    Rcpp::Named("scale_config") = s.config,
    Rcpp::Named("scale_target") = s.target,
    Rcpp::Named("n_actions_requested") = s.n_actions_requested,
    Rcpp::Named("n_actions_sampled") = s.n_actions_sampled,
    Rcpp::Named("n_moves_requested") = s.n_moves_requested,
    Rcpp::Named("n_moves_valid") = s.n_moves_valid,
    Rcpp::Named("n_cost_nonzero") = s.n_cost_nonzero,
    Rcpp::Named("n_config_nonzero") = s.n_config_nonzero,
    Rcpp::Named("n_target_nonzero") = s.n_target_nonzero,
    Rcpp::Named("calibrated") = s.calibrated,
    Rcpp::Named("by_action") = Rcpp::DataFrame::create(
      Rcpp::Named("action") = action,
      Rcpp::Named("n_eligible_cells") = n_eligible,
      Rcpp::Named("n_moves_requested") = n_requested,
      Rcpp::Named("n_moves_valid") = n_valid,
      Rcpp::Named("n_cost_nonzero") = n_cost_nonzero,
      Rcpp::Named("n_config_nonzero") = n_config_nonzero,
      Rcpp::Named("n_target_nonzero") = n_target_nonzero,
      Rcpp::Named("median_abs_cost") = median_abs_cost,
      Rcpp::Named("median_abs_config") = median_abs_config,
      Rcpp::Named("median_abs_target") = median_abs_target
    )
  );
}

} // namespace


//------------------------------- Main functions -------------------------------

// [[Rcpp::export]]
Rcpp::List run_ebrel_cpp(
    const std::vector<int>& E_int,
    const std::vector<double>& C,
    const std::vector<double>& SD,
    const std::vector<int>& D,
    const std::vector<double>& SxH,
    const std::vector<double>& O,
    const std::vector<double>& action_weight,
    const std::vector<int>& LM_int,
    const std::vector<int>& X0,
    int dim_x,
    int dim_y,
    int n_actions,
    int n_species,
    double sentinel,
    Rcpp::List opt
) {
  
  // ---- Number of cells ----
  const int n_cells = dim_x * dim_y;
  
  // ---- cast compact types ----
  std::vector<int8_t>  E(E_int.begin(), E_int.end());
  std::vector<uint8_t> LM(LM_int.begin(), LM_int.end());

  // ---- For build ebrel options ----
  const int    universal_disp_thres = opt_int(opt, "universal_disp_thres", 20);
  const int    max_disp_steps       = opt_int(opt, "max_disp_steps", 10);
  const int    roi_cap              = opt_int(opt, "roi_cap", 100);
  const bool   precompute_W         = opt_bool(opt, "precompute_W", true);
  
  // ---- Sigma_in: NA / missing => compute default in C++ ----
  const double sigma_in = opt_double(opt, "sigma_in", -1);

  // These are preference weights selected by the user. They are deliberately
  // not rescaled here.
  const double alpha = opt_double(opt, "alpha", 1.0);
  const double beta = opt_double(opt, "beta", 25.0);
  const double gamma = opt_double(opt, "gamma", 100.0);
  
  // ---- Convert optional supplied X0 ----
  std::vector<int8_t> X0_seed;
  
  if (!X0.empty()) {
    if (X0.size() != static_cast<std::size_t>(n_cells)) {
      throw std::runtime_error(
          "X0 must have length dim_x * dim_y."
      );
    }
    
    X0_seed.reserve(X0.size());
    
    for (const int value : X0) {
      if (value < -1 || value >= n_actions) {
        throw std::runtime_error(
            "X0 values must be -1 or valid zero-indexed action codes."
        );
      }
      
      X0_seed.push_back(static_cast<int8_t>(value));
    }
  }
  
  // ---- build fully-initialized input + defaults ----
  auto built = build_ebrel(
    std::move(E),
    std::vector<double>(C.begin(), C.end()),
    std::vector<double>(SD.begin(), SD.end()),
    std::vector<int>(D.begin(), D.end()),
    std::vector<double>(SxH.begin(), SxH.end()),
    std::vector<double>(O.begin(), O.end()),
    std::vector<float>(action_weight.begin(), action_weight.end()), 
    std::move(LM),
    std::move(X0_seed),
    dim_x, dim_y, 
    n_actions,
    n_species,
    sentinel,
    sigma_in,
    universal_disp_thres,
    max_disp_steps,
    roi_cap,
    precompute_W
  );
  
  // ---- Apply all options before any calibration/tuning that depends on them ----
  apply_run_options(opt, built.opt);
  const TempTuneOptions tune = parse_temp_tune_opts(opt);
  
  // ---- RunEBRELInput must contain these new fields after the refactor: ----
  //   double alpha, beta, gamma;
  //   ObjectiveScales objective_scales;
  built.in.alpha = alpha;
  built.in.beta = beta;
  built.in.gamma = gamma;
  
  // ---- Calculate improvement weighting for W ----
  const int improve_action = built.in.n_actions - 1; // zero-index
  
  int n_improvable = 0;
  
  for (int cell = 0; cell < n_cells; ++cell) {
    const int idx = action_cell_index(
      cell,
      improve_action,
      n_cells
    );
    
    if (built.in.U[static_cast<std::size_t>(idx)] == 0u) {
      ++n_improvable;
    }
  }
  
  built.in.improve_w_default = (n_improvable > 0)
    ? 1.0 / static_cast<double>(n_improvable)
    : 0.0;

  // ---- Generate X0 only when none was supplied ----
  if (built.in.X0.empty()) {
    const int x0_seed =
      (built.opt.rng_seed >= 0)
    ? built.opt.rng_seed
    : tune.rng_seed;
    
    built.in.X0 = generate_X0_CI(
      built.in.U,
      built.in.n_actions,
      built.in.dim_x,
      built.in.dim_y,
      built.opt.base_prob_X0,
      x0_seed
    );
  }
  
  // Whether supplied or generated, initialise derived X0 state here.
  built.in.improve_count = initialise_improve_count(
    built.in.X0,
    built.in.species_plan
  );
  
  // ---- Objective-scale calibration ----
  // This intentionally does not use X0. Each sample is an addition to an
  // empty planning solution, while compute_H retains the existing habitat /
  // species/dispersal context encoded in BuiltEBREL.
  const int n_scale_samples = opt_int(opt, "n_scale_samples", 500);
  const int calibration_seed = opt_int(
    opt,
    "calibration_seed",
    (built.opt.rng_seed >= 0) ? built.opt.rng_seed + 104729 : -1
  );
  std::mt19937 calibration_rng = make_rng(calibration_seed);
  
  const std::vector<int8_t> X_calibration(
      static_cast<std::size_t>(n_cells),
      static_cast<int8_t>(-1)
  );
  
  const RawObjectiveEvaluator evaluate_raw =
    [&](const std::vector<int8_t>& X) -> RawObjectiveComponents {
      const std::vector<int> improve_count = initialise_improve_count(
        X,
        built.in.species_plan
      );
      
      const HResult h = compute_H(
        X,
        built.in.C,
        1.0, 1.0, 1.0,
        1.0, 1.0, 1.0,
        built.in.n_actions,
        built.in.n_species,
        built.in.n_habitats,
        built.in.dim_x,
        built.in.dim_y,
        built.in.universal_disp_thres,
        built.in.max_disp_steps,
        built.in.roi_cap,
        built.in.LM,
        built.in.row_first_land,
        built.in.row_last_land,
        built.in.col_first_land,
        built.in.col_last_land,
        built.in.E,
        built.in.Etiles_per_h,
        built.in.cell_r,
        built.in.cell_c,
        built.in.rowruns_cache,
        built.in.species_plan,
        improve_count
      );
      
      static int calibration_debug_calls = 0;
      
      if (calibration_debug_calls < 20) {
        
        int selected_cell = -1;
        int selected_action = -1;
        
        for (int cell = 0; cell < n_cells; ++cell) {
          if (X[static_cast<std::size_t>(cell)] >= 0) {
            selected_cell = cell;
            selected_action = X[static_cast<std::size_t>(cell)];
            break;
          }
        }
        
        std::cout
        << "[calibration evaluator]"
        << " call=" << calibration_debug_calls
        << " selected_cell=" << selected_cell
        << " selected_action=" << selected_action
        << " h.F1=" << h.F1
        << " h.F2=" << h.F2
        << " h.gx=" << h.gx
        << '\n';
        
        for (int sp = 0; sp < built.in.n_species; ++sp) {
          const std::size_t spi = static_cast<std::size_t>(sp);
          
          std::cout
          << "  sp=" << sp
          << " active="
          << static_cast<int>(built.in.species_plan.active[spi])
          << " seed_len=" << built.in.species_plan.seed_len[spi]
          << " O_n=" << built.in.species_plan.O_n[spi]
          << " improve_count=" << improve_count[spi]
          << " g=" << h.g[spi]
          << " g_create=" << h.g_create[spi]
          << " g_improve=" << h.g_improve[spi]
          << '\n';
        }
      }
      
      ++calibration_debug_calls;
      
      return RawObjectiveComponents{
        h.F1,
        h.F2,
        h.gx
      };
    };
    
    const ObjectiveScalingCalibration calibration =
      calibrate_objective_scales(
        built.in.U,
        X_calibration,
        built.in.n_actions,
        built.in.dim_x,
        built.in.dim_y,
        n_scale_samples,
        evaluate_raw,
        calibration_rng
      );
    
    built.in.objective_scales = calibration.scales;
    
    if (built.opt.verbose) {
      std::cout
      << "[objective calibration]"
      << " scale_cost=" << calibration.scales.cost
      << " scale_config=" << calibration.scales.config
      << " scale_target=" << calibration.scales.target
      << " valid_moves=" << calibration.scales.n_moves_valid
      << '\n';
      
      std::cout
      << "[objective preferences]"
      << " alpha=" << built.in.alpha
      << " beta=" << built.in.beta
      << " gamma=" << built.in.gamma
      << '\n';
    }
  
    // ---- Optional tuning of initial temperature ----
    Rcpp::List tune_diag = Rcpp::List::create(
      Rcpp::Named("enabled") = false
    );
    
    if (tune.enabled) {
      // The temperature tuner must accept user weights and fixed objective
      // scales, and combine them exactly as the SA loop does.
      const InitTempResult t0 = estimate_initial_temperature_benameur_cpp(
        built.in.X0,
        built.in.W,
        built.in.U,
        built.in.C,
        built.in.E,
        built.in.Etiles_per_h,
        built.in.cell_r,
        built.in.cell_c,
        built.in.rowruns_cache,
        built.in.species_plan,
        built.in.universal_disp_thres,
        built.in.max_disp_steps,
        built.in.roi_cap,
        built.in.LM,
        built.in.row_first_land,
        built.in.row_last_land,
        built.in.col_first_land,
        built.in.col_last_land,
        built.in.n_actions,
        built.in.n_species,
        built.in.n_habitats,
        built.in.dim_x,
        built.in.dim_y,
        built.in.alpha,
        built.in.beta,
        built.in.gamma,
        built.in.objective_scales.cost,
        built.in.objective_scales.config,
        built.in.objective_scales.target,
        built.opt.step_proportion,
        built.opt.step_probability,
        tune.num_samples,
        tune.chi0,
        tune.p,
        tune.tol_logchi,
        tune.max_iters,
        tune.T1,
        tune.max_tries_factor,
        tune.rng_seed,
        tune.verbose
      );
      
      built.opt.temp = t0.T0;
      
      tune_diag = Rcpp::List::create(
        Rcpp::Named("enabled") = true,
        Rcpp::Named("T0") = t0.T0,
        Rcpp::Named("chi_hat_final") = t0.diag.chi_hat_final,
        Rcpp::Named("chi0") = t0.diag.chi0,
        Rcpp::Named("iters") = t0.diag.iters,
        Rcpp::Named("samples_used") = t0.diag.samples_used,
        Rcpp::Named("tries") = t0.diag.tries,
        Rcpp::Named("dE_median") = t0.diag.dE_median,
        Rcpp::Named("dE_mean") = t0.diag.dE_mean
      );
    }
    
    // ---- Run EBREL -------------------------------------------------------------
    const RunEBRELResult res = run_ebrel(built.in, built.opt);
    
    Rcpp::IntegerVector X_best_iv(res.X_best.size());
    for (std::size_t i = 0; i < res.X_best.size(); ++i) {
      X_best_iv[i] = static_cast<int>(res.X_best[i]);
    }
    
    Rcpp::IntegerVector X0_iv(built.in.X0.size());
    for (std::size_t i = 0; i < built.in.X0.size(); ++i) {
      X0_iv[i] = static_cast<int>(built.in.X0[i]);
    }
    
    Rcpp::IntegerVector U_iv(built.in.U.size());
    for (std::size_t i = 0; i < built.in.U.size(); ++i) {
      U_iv[i] = static_cast<int>(built.in.U[i]);
    }
    U_iv.attr("dim") = Rcpp::IntegerVector::create(
      n_cells,
      built.in.n_actions
    );
    
    return Rcpp::List::create(
      Rcpp::Named("X0") = X0_iv,
      Rcpp::Named("U") = U_iv,
      Rcpp::Named("X_best") = X_best_iv,
      Rcpp::Named("H_best") = res.H_best,
      Rcpp::Named("iterations_run") = res.iterations_run,
      Rcpp::Named("accepted") = res.accepted,
      Rcpp::Named("proposals") = res.proposals,
      Rcpp::Named("overall_acc") = res.overall_acc,
      Rcpp::Named("early_stop_iter") = res.early_stop_iter,
      Rcpp::Named("H_trace") = res.H_trace,
      Rcpp::Named("F_trace") = res.F_trace,
      Rcpp::Named("F1_trace") = res.F1_trace,
      Rcpp::Named("F2_trace") = res.F2_trace,
      Rcpp::Named("g_best") = res.g_best,
      Rcpp::Named("g_create_best") = res.g_create_best,
      Rcpp::Named("g_improve_best") = res.g_improve_best,
      Rcpp::Named("acc_rate_trace") = res.acc_rate_trace,
      Rcpp::Named("iter_ms_total") = res.iter_ms_total,
      Rcpp::Named("iter_count") = res.iter_count,
      Rcpp::Named("objective_weights") = Rcpp::List::create(
        Rcpp::Named("alpha") = built.in.alpha,
        Rcpp::Named("beta") = built.in.beta,
        Rcpp::Named("gamma") = built.in.gamma
      ),
      Rcpp::Named("objective_calibration") = as_r_calibration(calibration),
      Rcpp::Named("temperature_tuning") = tune_diag
    );
}


// [[Rcpp::export]]
Rcpp::IntegerVector generate_X0_CI_R(const std::vector<uint8_t>& U,
                                    int n_actions,
                                    int dim_x,
                                    int dim_y,
                                    double base_prob,
                                    int seed) {
  std::vector<int8_t> x0 = generate_X0_CI(U, n_actions, dim_x, dim_y, base_prob, seed);

  Rcpp::IntegerVector out(x0.size());
  for (std::size_t i = 0; i < x0.size(); ++i) out[i] = static_cast<int>(x0[i]);

  return out;
}
