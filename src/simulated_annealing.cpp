//------------------------------ Simulated annealing ---------------------------

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <chrono>
#include <fstream>

#include "simulated_annealing.h"
#include "dispersal_utils.h"
#include "update_candidate.h"
#include "objective_utils.h"
#include "optimisation_utils.h"
#include "species_plan.h"

/* Simulated annealing algorithm with:
 - option for Lam et al style tuning schedule [Lam, J. and Delosme, J.M., 1988.
 An efficient simulated annealing schedule: derivation. New Haven, CT: Yale
 Electrical Engineering Department, 8816.]
 - or harmonic/geometric cooling schedule with a option to control rate of cooling
 */
SAResult simulated_annealing(
    // --------- Input data and parameters
    const std::vector<int8_t>& X0,
    const std::vector<double>& W,
    const std::vector<uint8_t>& U,
    const std::vector<double>& C,
    const std::vector<int8_t>& E,
    const std::vector<std::vector<std::size_t>>& Etiles_per_h,
    const std::vector<int>&    cell_r,
    const std::vector<int>&    cell_c,
    const RowRunsCache& rowruns_cache,
    const SpeciesPlan& species_plan,
    std::vector<int> improve_count,
    // ---------- Dimensions
    int n_actions,
    int n_species,
    int n_habitats,
    int dim_x,
    int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    // ---------- Cap for size of ROI (in grid cells)
    int roi_cap,
    // ---------- For avoiding extra work on sea cells
    const std::vector<uint8_t>& LM,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
    // ---------- Objective function related
    double alpha,
    double beta,
    double gamma,
    // ---------- Internal raw-component standardisation constants
    double scale_cost,
    double scale_config,
    double scale_target,
    // ---------- Scaling for improvement sampling
    double improve_total_weight,
    // ---------- Simulated Annealing
    double step_proportion,
    double step_probability,
    int n_iterations,
    double temp,
    double cooling_rate_c,        // constant for tuning cooling rate
    // ---------- LAM: user controls (default = disabled) ----------
    bool   lam_enabled,           // turn Lam-style online control on/off
    double lam_target_mid,        // target uphill acceptance during early/mid run
    double lam_target_final,      // target uphill acceptance at the end
    double lam_hold_frac,         // fraction of run to hold lam_target_mid before decaying
    double lam_p,                 // damping exponent in Ben-Ameur correction (>=1)
    // ---------- Early stopping
    int min_iterations,           // require at least this many iterations
    int acceptance_window,        // window length for acceptance rate
    double acceptance_thres,      // "low" acceptance threshold
    int iter_no_improve,          // consecutive iters with no meaningful improvement
    double improve_eps,           // relative improvement needed to reset patience
    // --------- Output controls
    int  write_every,             // 0 = disabled
    std::string trace_file,       // empty = disabled
    bool verbose
) {

  // Printing control
  std::ios_base::sync_with_stdio(false);
  std::cout.tie(nullptr);

  // Write out
  std::ofstream trace_out;
  if (write_every > 0 && !trace_file.empty()) {
    trace_out.open(trace_file, std::ios::out | std::ios::trunc);
    if (!trace_out) {
      throw std::runtime_error("Failed to open trace file: " + trace_file);
    }

    // Header (CSV-style)
    trace_out << "iter,best_H,curr_H,temp,accepted_total,attempted_total,overall_acc,avg_ms\n";
    trace_out.flush();
  }

  // RNG
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);

  // Histories
  std::vector<double> H_history;
  std::vector<double> F_history;
  std::vector<double> F1_history;
  std::vector<double> F2_history;

  // Initial evaluation of X0
  HResult init_scores = compute_H(
    X0, C,
    alpha, beta, gamma,
    scale_cost, scale_config, scale_target,
    n_actions, n_species, n_habitats,
    dim_x, dim_y,
    universal_disp_thres, max_disp_steps, roi_cap,
    LM,
    row_first_land, row_last_land,
    col_first_land, col_last_land,
    E,
    Etiles_per_h,
    cell_r, cell_c,
    rowruns_cache,
    species_plan,
    improve_count
  );

  double best_score = init_scores.H;
  double best_Fx    = init_scores.Fx;
  //double best_F1    = init_scores.F1;
  //double best_F2    = init_scores.F2;
  
  // Current & best states
  std::vector<int8_t> best = X0;
  std::vector<int8_t> curr = X0;
  double curr_eval = best_score;

  std::vector<double> g_best = init_scores.g;
  std::vector<double> g_create_best = init_scores.g_create;
  std::vector<double> g_improve_best = init_scores.g_improve;

  // Ensure stopping parameters are in a sensible range
  const int    min_iters    = std::max(0, min_iterations);
  const int    win          = std::max(1, acceptance_window);
  const double acc_min      = std::clamp(acceptance_thres, 0.0, 1.0);
  const int    patience     = std::max(1, iter_no_improve);
  const double improve_thr  = std::max(0.0, improve_eps);

  // Diagnostics
  std::vector<double> acc_rate_trace; // per window
  int    early_stop_iter = -1;
  int    attempted_total = 0;
  int    accepted_total  = 0;

  // For stopping criterion
  int    attempted_in_win = 0;
  int    accepted_in_win  = 0;
  int    no_improve       = 0;
  double last_best        = best_score;

  // Timings
  long long iter_ms_total = 0;
  int       iter_count    = 0;

  // ---------- LAM: helpers & counters ----------
  auto clamp01 = [](double x, double eps = 1e-12) {
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
  };
  auto target_acceptance = [&](double progress01) {
    progress01 = std::min(1.0, std::max(0.0, progress01));
    if (progress01 <= lam_hold_frac) return lam_target_mid;
    double u = (progress01 - lam_hold_frac) / (1.0 - lam_hold_frac);
    return std::exp((1.0 - u) * std::log(lam_target_mid) + u * std::log(lam_target_final));
  };
  int uphill_attempted_in_win = 0;
  int uphill_accepted_in_win  = 0;
  double temp_curr = temp; // for LAM

  if (verbose) {
    std::cout
      << "  Fx=" << init_scores.Fx
      << ", F1=" << init_scores.F1
      << ", F2=" << init_scores.F2
      << ", alpha*(F1/scale_cost)=" << alpha * (init_scores.F1 / scale_cost)
      << ", beta*(F2/scale_config)=" << beta * (init_scores.F2 / scale_config)
      << ", gamma*(gx/scale_target)=" << gamma * (init_scores.gx / scale_target)
      << ", H=" << init_scores.H
      << '\n';
  }
  
  // Main loop
  for (int z = 0; z < n_iterations; ++z) {

    // Propose candidate
    const uint32_t rng_seed = rng(); // pulls from SA RNG
    CandidateUpdate update = update_candidate(
      W, U, curr,
      step_proportion,
      step_probability,
      improve_total_weight,
      n_actions, 
      n_habitats,
      dim_x, dim_y,
      rng_seed
    );

    if (update.changes.empty()) continue;
    const std::vector<int8_t>& candidate = update.X;

    // Tracking acceptance
    attempted_in_win++;
    attempted_total++;
  
    // Save improve_count before proposal
    std::vector<int> improve_count_old = improve_count;
    std::vector<int> improve_count_candidate = improve_count;

    const int improve_action = species_plan.n_habitats;

    for (const auto& ch : update.changes) {
      const bool old_improve = (static_cast<int>(ch.old_action) == improve_action);
      const bool new_improve = (static_cast<int>(ch.new_action) == improve_action);
      
      if (old_improve == new_improve) continue;
      
      const int delta = new_improve ? 1 : -1;
      const uint32_t start = species_plan.cell_species_start[ch.cell];
      const uint32_t len   = species_plan.cell_species_len[ch.cell];
      
      for (uint32_t j = 0; j < len; ++j) {
        improve_count_candidate[species_plan.cell_species_pool[start + j]] += delta;
      }
    }

    // Evaluate candidate with new candidate
    auto t0 = std::chrono::steady_clock::now();
    HResult scores = compute_H(
      candidate, C,
      alpha, beta, gamma,
      scale_cost, scale_config, scale_target,
      n_actions, n_species, n_habitats,
      dim_x, dim_y,
      universal_disp_thres, max_disp_steps, roi_cap,
      LM,
      row_first_land, row_last_land,
      col_first_land, col_last_land,
      E,
      Etiles_per_h,
      cell_r, cell_c,
      rowruns_cache,
      species_plan,
      improve_count_candidate
    );
    double dt = ms_since(t0);
    iter_ms_total += static_cast<long long>(dt);
    iter_count    += 1;
    
    // Candidate evaluation
    double candidate_eval = scores.H;

    // Trace histories
    H_history.push_back(candidate_eval);
    F_history.push_back(scores.Fx);
    F1_history.push_back(scores.F1);
    F2_history.push_back(scores.F2);

    // Update best
    if (candidate_eval < best_score) {
      best = candidate;
      best_score = candidate_eval;
      best_Fx = scores.Fx;
      //best_F1 = scores.F1;
      //best_F2 = scores.F2;
      g_best = scores.g;
      g_create_best = scores.g_create;
      g_improve_best = scores.g_improve;
    }

    // ------------- Metropolis acceptance -------------------------------
    const double delta = candidate_eval - curr_eval;
    const double epsE  = 1e-12 * (1.0 + std::max(std::abs(curr_eval), std::abs(candidate_eval)));

    const bool downhill = (delta < -epsE);
    const bool uphill   = (delta >  epsE);
    const bool flat     = !downhill && !uphill;

    if (uphill) uphill_attempted_in_win++;

    const double t = (lam_enabled)
      ? temp_curr
      : (temp / (1.0 + cooling_rate_c * z));  // harmonic default

    bool accepted = false;
    if (downhill || flat) {
      accepted = true;
    } else if (t > 0.0 && std::isfinite(t)) {
      double u = runif(rng);
      if (u <= 0.0) u = std::numeric_limits<double>::min();
      const double logu = std::log(u);
      const double logA = -delta / t;
      accepted = (logu < std::min(0.0, logA));
    }

    if (accepted) {
      curr = candidate;
      curr_eval = candidate_eval;
      improve_count = std::move(improve_count_candidate);
      accepted_in_win++;
      accepted_total++;
      if (uphill) uphill_accepted_in_win++;
    } else {
      improve_count = improve_count_old;
    }

    // ---- Stopping logic ----
    {
      double denom = std::max(1.0, std::abs(last_best));
      double rel_delta = std::abs(best_score - last_best) / denom;
      if (rel_delta >= improve_thr) {
        no_improve = 0;
        last_best = best_score;
      } else {
        no_improve++;
      }
    }

    if ((z + 1) % win == 0) {

      double acc_rate = (attempted_in_win > 0)
      ? static_cast<double>(accepted_in_win) / attempted_in_win
      : 0.0;

      acc_rate_trace.push_back(acc_rate);

      if ((z + 1) >= min_iters && acc_rate < acc_min && no_improve >= patience) {
        std::cout << "Stopping early at iteration " << (z + 1)
                  << " due to combined criterion: low acceptance AND no improvement.\n";
        early_stop_iter = (z + 1);
        break;
      }

      if (lam_enabled) {
        double hat_up = (uphill_attempted_in_win > 0)
        ? static_cast<double>(uphill_accepted_in_win) / uphill_attempted_in_win
        : ((attempted_in_win > 0) ? static_cast<double>(accepted_in_win) / attempted_in_win : 0.0);

        double progress = static_cast<double>(z + 1) / static_cast<double>(n_iterations);
        double a_star   = target_acceptance(progress);

        double num = std::log(clamp01(hat_up));
        double den = std::log(clamp01(a_star));
        double ratio = (den != 0.0) ? (num / den) : 1.0;
        if (!(ratio > 0.0) || !std::isfinite(ratio)) ratio = 1.0;

        temp_curr *= std::pow(ratio, 1.0 / lam_p);
        if (!std::isfinite(temp_curr) || temp_curr < 1e-9) temp_curr = 1e-9;
        if (temp_curr > 1e16) temp_curr = 1e16;
      }

      // reset window counters
      attempted_in_win = 0;
      accepted_in_win  = 0;
      uphill_attempted_in_win = 0;
      uphill_accepted_in_win  = 0;
    }

    // -------- Write intermediate results to file --------
    if (write_every > 0 && trace_out.is_open() &&
        ((z + 1) % write_every == 0)) {

      double overall_acc =
        (attempted_total > 0)
      ? static_cast<double>(accepted_total) / attempted_total
      : std::numeric_limits<double>::quiet_NaN();

      double T_curr = lam_enabled
      ? temp_curr
      : (temp / (1.0 + cooling_rate_c * z));

      double avg_ms =
        (iter_count > 0)
        ? static_cast<double>(iter_ms_total) / iter_count
      : 0.0;

      trace_out << (z + 1) << ","
                << best_score << ","
                << curr_eval << ","
                << T_curr << ","
                << accepted_total << ","
                << attempted_total << ","
                << overall_acc << ","
                << avg_ms
                << '\n';

      // Optional but recommended for long runs
      trace_out.flush();
    }

    if (verbose && ((z + 1) % win == 0)) {
      std::cout << "[SA] it=" << (z + 1)
                << " T=" << (lam_enabled ? temp_curr : (temp / (1.0 + cooling_rate_c * z)))
                << " Best H(x) = " << best_score
                << " Best F(x) = " << best_Fx
                << '\n' << std::flush;
    }
  }

  // Result of SA
  SAResult out;
  out.X_best   = best;
  out.H_best   = best_score;
  out.H_trace  = H_history;
  out.F_trace  = F_history;
  out.F1_trace = F1_history;
  out.F2_trace = F2_history;
  out.g_best   = g_best;
  out.g_create_best   = g_create_best;
  out.g_improve_best   = g_improve_best;

  // Diagnostics
  out.diag.acceptance_rate_trace = std::move(acc_rate_trace);
  out.diag.early_stop_iter       = early_stop_iter;
  out.diag.attempted_total       = attempted_total;
  out.diag.accepted_total        = accepted_total;
  out.diag.iter_ms_total         = iter_ms_total;
  out.diag.iter_count            = iter_count;

  const int proposals_from_trace = static_cast<int>(H_history.size());
  const int proposals_print = (attempted_total > 0) ? attempted_total : proposals_from_trace;

  double overall_acc = (proposals_print > 0)
    ? static_cast<double>(accepted_total) / proposals_print
  : std::numeric_limits<double>::quiet_NaN();

  std::cout << "[SA] proposals=" << proposals_print
            << " | accepted=" << accepted_total
            << " | overall_acc=" << overall_acc
            << (early_stop_iter > 0 ? " | early_stop_at=" + std::to_string(early_stop_iter) : "")
            << std::endl;

  return out;
}
