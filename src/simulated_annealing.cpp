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
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "simulated_annealing.h"
#include "dispersal_utils.h"
#include "update_candidate.h"
#include "objective_utils.h"
#include "optimisation_utils.h"
#include "species_plan.h"

//-------------------------- Local helpers -------------------------------------

namespace {

template <typename T>
void write_binary(
    std::ofstream& out,
    const T& value
) {
  out.write(
    reinterpret_cast<const char*>(&value),
    sizeof(T)
  );
  
  if (!out) {
    throw std::runtime_error(
        "Failed while writing checkpoint"
    );
  }
}

template <typename T>
void read_binary(
    std::ifstream& in,
    T& value
) {
  in.read(
    reinterpret_cast<char*>(&value),
    sizeof(T)
  );
  
  if (!in) {
    throw std::runtime_error(
        "Failed while reading checkpoint"
    );
  }
}


template <typename T>
void write_vector(
    std::ofstream& out,
    const std::vector<T>& x
) {
  const std::uint64_t n =
    static_cast<std::uint64_t>(x.size());
  
  write_binary(out, n);
  
  if (n > 0) {
    out.write(
      reinterpret_cast<const char*>(x.data()),
      static_cast<std::streamsize>(
        n * sizeof(T)
      )
    );
    
    if (!out) {
      throw std::runtime_error(
          "Failed while writing checkpoint vector"
      );
    }
  }
}


template <typename T>
void read_vector(
    std::ifstream& in,
    std::vector<T>& x
) {
  std::uint64_t n = 0;
  
  read_binary(in, n);
  
  x.resize(
    static_cast<std::size_t>(n)
  );
  
  if (n > 0) {
    in.read(
      reinterpret_cast<char*>(x.data()),
      static_cast<std::streamsize>(
        n * sizeof(T)
      )
    );
    
    if (!in) {
      throw std::runtime_error(
          "Failed while reading checkpoint vector"
      );
    }
  }
}


void write_string(
    std::ofstream& out,
    const std::string& x
) {
  const std::uint64_t n =
    static_cast<std::uint64_t>(x.size());
  
  write_binary(out, n);
  
  if (n > 0) {
    out.write(
      x.data(),
      static_cast<std::streamsize>(n)
    );
    
    if (!out) {
      throw std::runtime_error(
          "Failed while writing checkpoint string"
      );
    }
  }
}


void read_string(
    std::ifstream& in,
    std::string& x
) {
  std::uint64_t n = 0;
  
  read_binary(in, n);
  
  x.resize(
    static_cast<std::size_t>(n)
  );
  
  if (n > 0) {
    in.read(
      x.data(),
      static_cast<std::streamsize>(n)
    );
    
    if (!in) {
      throw std::runtime_error(
          "Failed while reading checkpoint string"
      );
    }
  }
}

}

// Write sa checkpoint data for restart
void write_sa_checkpoint(
    const std::string& checkpoint_file,
    const SACheckpoint& cp
) {
  
  const std::string tmp_file =
    checkpoint_file + ".tmp";
  
  std::ofstream out(
      tmp_file,
      std::ios::binary |
        std::ios::trunc
  );
  
  if (!out) {
    throw std::runtime_error(
        "Failed to open checkpoint file for writing: " +
          tmp_file
    );
  }
  
  // File identity / version
  
  const std::uint32_t magic =
    0x45425232; // "EBR2"
  
  const std::uint32_t version = 2;
  
  write_binary(out, magic);
  write_binary(out, version);

  // Metadata
  write_binary(out, cp.dim_x);
  write_binary(out, cp.dim_y);
  write_binary(out, cp.n_actions);
  write_binary(out, cp.n_species);
  write_binary(out, cp.n_habitats);
  write_binary(out, cp.n_iterations);
  
  // Objective scaling
  write_binary(out, cp.scale_cost);
  write_binary(out, cp.scale_config);
  write_binary(out, cp.scale_target);
  
  // Main state
  write_binary(out, cp.iteration);
  write_vector(out, cp.curr);
  write_binary(out, cp.curr_eval);
  write_vector(out, cp.improve_count);
  write_vector(out, cp.best);
  write_binary(out, cp.best_score);
  write_binary(out, cp.best_Fx);
  write_vector(out, cp.g_best);
  write_vector(out, cp.g_create_best);
  write_vector(out, cp.g_improve_best);

  // Counters
  write_binary(out, cp.attempted_total);
  write_binary(out, cp.accepted_total);
  write_binary(out, cp.attempted_in_win);
  write_binary(out, cp.accepted_in_win);
  write_binary(out, cp.no_improve);
  write_binary(out, cp.last_best);
  write_binary(out, cp.uphill_attempted_in_win);
  write_binary(out, cp.uphill_accepted_in_win);
  write_binary(out, cp.temp_curr);

  // Timing
  write_binary(out, cp.iter_ms_total);
  write_binary(out, cp.iter_count);
  
  // Histories
  write_vector(out, cp.H_history);
  write_vector(out, cp.F_history);
  write_vector(out, cp.F1_history);
  write_vector(out, cp.F2_history);
  write_vector(out, cp.acc_rate_trace);
  
  // RNG
  write_string(out, cp.rng_state);
  
  out.flush();
  out.close();
  
  if (!out) {
    throw std::runtime_error(
        "Failed to complete checkpoint write: " +
          tmp_file
    );
  }
  
  // Atomic-ish replacement
  std::error_code ec;
  std::filesystem::rename(
    tmp_file,
    checkpoint_file,
    ec
  );
  
  if (ec) {
    
    // Destination may already exist
    std::filesystem::remove(
      checkpoint_file,
      ec
    );
    
    ec.clear();
    
    std::filesystem::rename(
      tmp_file,
      checkpoint_file,
      ec
    );
  }
  
  if (ec) {
    throw std::runtime_error(
        "Failed to move checkpoint into place: " +
          checkpoint_file +
          " (" + ec.message() + ")"
    );
  }
}

// Read sa checkpoint data for restart
SACheckpoint read_sa_checkpoint(
    const std::string& checkpoint_file
) {
  
  std::ifstream in(
      checkpoint_file,
      std::ios::binary
  );
  
  if (!in) {
    throw std::runtime_error(
        "Failed to open checkpoint file: " +
          checkpoint_file
    );
  }
  
  SACheckpoint cp;
  
  
  // Validate file identity
  
  std::uint32_t magic = 0;
  std::uint32_t version = 0;
  
  read_binary(in, magic);
  read_binary(in, version);
  
  constexpr std::uint32_t expected_magic =
    0x45425232;
  
  if (magic != expected_magic) {
    throw std::runtime_error(
        "File is not a valid ebrel2 SA checkpoint: " +
          checkpoint_file
    );
  }
  
  if (version != 2) {
    throw std::runtime_error(
        "Unsupported SA checkpoint version: " +
          std::to_string(version)
    );
  }
  
  
  // Metadata
  read_binary(in, cp.dim_x);
  read_binary(in, cp.dim_y);
  read_binary(in, cp.n_actions);
  read_binary(in, cp.n_species);
  read_binary(in, cp.n_habitats);
  read_binary(in, cp.n_iterations);
  
  // Objective scaling
  read_binary(in, cp.scale_cost);
  read_binary(in, cp.scale_config);
  read_binary(in, cp.scale_target);
  
  // Main state
  read_binary(in, cp.iteration);
  read_vector(in, cp.curr);
  read_binary(in, cp.curr_eval);
  read_vector(in, cp.improve_count);
  read_vector(in, cp.best);
  read_binary(in, cp.best_score);
  read_binary(in, cp.best_Fx);
  read_vector(in, cp.g_best);
  read_vector(in, cp.g_create_best);
  read_vector(in, cp.g_improve_best);
  
  // Counters
  read_binary(in, cp.attempted_total);
  read_binary(in, cp.accepted_total);
  read_binary(in, cp.attempted_in_win);
  read_binary(in, cp.accepted_in_win);
  read_binary(in, cp.no_improve);
  read_binary(in, cp.last_best);
  read_binary(in, cp.uphill_attempted_in_win);
  read_binary(in, cp.uphill_accepted_in_win);
  read_binary(in, cp.temp_curr);

  // Timing
  read_binary(in, cp.iter_ms_total);
  read_binary(in, cp.iter_count);
    
  // Histories
  read_vector(in, cp.H_history);
  read_vector(in, cp.F_history);
  read_vector(in, cp.F1_history);
  read_vector(in, cp.F2_history);
  read_vector(in, cp.acc_rate_trace);

  // RNG
  read_string(in, cp.rng_state);
  
  if (!in) {
    throw std::runtime_error(
        "Checkpoint file is incomplete or corrupt: " +
          checkpoint_file
    );
  }
  
  return cp;
}

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
    // --------- Output / restart controls
    int trace_every,                        // 0 = disabled
    const std::string& trace_file,          // empty = disabled
    int checkpoint_every,                   // 0 = disabled
    const std::string& checkpoint_file,     // where checkpoints are written
    const std::string& restart_file,        // empty = new run
    bool verbose
) {

  // Printing control
  std::ios_base::sync_with_stdio(false);
  std::cout.tie(nullptr);

  // Write out
  std::ofstream trace_out;
  
  if (trace_every > 0 && !trace_file.empty()) {
    
    if (restart_file.empty()) {
      
      // New run: create a fresh trace
      trace_out.open(
        trace_file,
        std::ios::out | std::ios::trunc
      );
      
    } else {
      
      // Restart: append to existing trace
      trace_out.open(
        trace_file,
        std::ios::out | std::ios::app
      );
    }
    
    if (!trace_out) {
      throw std::runtime_error(
          "Failed to open trace file: " +
            trace_file
      );
    }
    
    // Only write the header for a new run
    if (restart_file.empty()) {
      trace_out
      << "iter,best_H,curr_H,temp,"
      << "accepted_total,attempted_total,"
      << "overall_acc,avg_ms\n";
      
      trace_out.flush();
    }
  }

  // Initialise iteration sequence
  int start_iter = 0;
  
  // RNG
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);

  // Histories
  std::vector<double> H_history;
  std::vector<double> F_history;
  std::vector<double> F1_history;
  std::vector<double> F2_history;
  std::vector<double> acc_rate_trace;
  
  // Mutable SA state: initialised either from X0 or checkpoint
  std::vector<int8_t> best;
  std::vector<int8_t> curr;
  
  double best_score;
  double best_Fx;
  double curr_eval;
  
  std::vector<double> g_best;
  std::vector<double> g_create_best;
  std::vector<double> g_improve_best;
  
  // Diagnostics
  int early_stop_iter = -1;
  int attempted_total = 0;
  int accepted_total  = 0;
  
  // For stopping criterion
  int attempted_in_win = 0;
  int accepted_in_win  = 0;
  int no_improve = 0;
  double last_best;
  
  // Timings
  long long iter_ms_total = 0;
  int iter_count = 0;
  
  int uphill_attempted_in_win = 0;
  int uphill_accepted_in_win  = 0;
  double temp_curr = temp;
  
  
  if (restart_file.empty()) {
    
    // =========================================
    // NEW RUN
    // =========================================
    
    HResult init_scores = compute_H(
      X0, C,
      alpha, beta, gamma,
      scale_cost, scale_config, scale_target,
      n_actions, n_species, n_habitats,
      dim_x, dim_y,
      universal_disp_thres,
      max_disp_steps,
      roi_cap,
      LM,
      row_first_land,
      row_last_land,
      col_first_land,
      col_last_land,
      E,
      Etiles_per_h,
      cell_r,
      cell_c,
      rowruns_cache,
      species_plan,
      improve_count
    );
    
    best_score = init_scores.H;
    best_Fx    = init_scores.Fx;
    
    best = X0;
    curr = X0;
    
    curr_eval = best_score;
    
    g_best         = init_scores.g;
    g_create_best  = init_scores.g_create;
    g_improve_best = init_scores.g_improve;
    
    last_best = best_score;
    
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
    
  } else {
    
    // =========================================
    // RESTART
    // =========================================
    
    // Load checkpoint data
    SACheckpoint cp =
      read_sa_checkpoint(restart_file);
    
  
    if (
        cp.dim_x != dim_x ||
          cp.dim_y != dim_y ||
          cp.n_actions != n_actions ||
          cp.n_species != n_species ||
          cp.n_habitats != n_habitats
    ) {
      throw std::runtime_error(
          "Checkpoint does not match current optimisation problem"
      );
    }
    
    // Restore objective scaling used by the original run
    scale_cost   = cp.scale_cost;
    scale_config = cp.scale_config;
    scale_target = cp.scale_target;
    
    // Restart iterations
    start_iter = cp.iteration;
    
    // Evaluations
    curr          = std::move(cp.curr);
    curr_eval     = cp.curr_eval;
    improve_count = std::move(cp.improve_count);
    
    best       = std::move(cp.best);
    best_score = cp.best_score;
    best_Fx    = cp.best_Fx;
    
    g_best         = std::move(cp.g_best);
    g_create_best  = std::move(cp.g_create_best);
    g_improve_best = std::move(cp.g_improve_best);
    
    attempted_total = cp.attempted_total;
    accepted_total  = cp.accepted_total;
    
    attempted_in_win = cp.attempted_in_win;
    accepted_in_win  = cp.accepted_in_win;
    
    no_improve = cp.no_improve;
    last_best  = cp.last_best;
    
    uphill_attempted_in_win =
      cp.uphill_attempted_in_win;
    
    uphill_accepted_in_win =
      cp.uphill_accepted_in_win;
    
    temp_curr = cp.temp_curr; // for LAM
    
    iter_ms_total = cp.iter_ms_total;
    iter_count    = cp.iter_count;
    
    H_history  = std::move(cp.H_history);
    F_history  = std::move(cp.F_history);
    F1_history = std::move(cp.F1_history);
    F2_history = std::move(cp.F2_history);
    
    acc_rate_trace =
      std::move(cp.acc_rate_trace);
    
    std::istringstream rng_in(cp.rng_state);
    
    rng_in >> rng >> runif;
    
    if (!rng_in) {
      throw std::runtime_error(
          "Failed to restore RNG state"
      );
    }
    
    if (verbose) {
      std::cout
      << "[SA] restarting from iteration "
      << start_iter
      << " | H_best=" << best_score
      << " | H_current=" << curr_eval
      << '\n';
    }
  }

  // Ensure stopping parameters are in a sensible range
  const int    min_iters    = std::max(0, min_iterations);
  const int    win          = std::max(1, acceptance_window);
  const double acc_min      = std::clamp(acceptance_thres, 0.0, 1.0);
  const int    patience     = std::max(1, iter_no_improve);
  const double improve_thr  = std::max(0.0, improve_eps);

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

  // Iteration completed counter
  int iterations_completed = start_iter;

  // Main loop
  for (int z = start_iter; z < n_iterations; ++z) {
    
    bool stop_now = false;
    
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

    // Only evaluate if a genuine proposal was generated
    if (!update.changes.empty()) {
      
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
        universal_disp_thres, 
        max_disp_steps, 
        roi_cap,
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

      if (uphill) {
        uphill_attempted_in_win++;
      }
      
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
      
      // Update iterations counter
      iterations_completed = z + 1;
    }

    // Acceptance window / Lam update
    if ((z + 1) % win == 0) {

      const double acc_rate =
        (attempted_in_win > 0)
      ? static_cast<double>(
          accepted_in_win
      ) / attempted_in_win
      : 0.0;

      acc_rate_trace.push_back(acc_rate);

      if (
          (z + 1) >= min_iters && 
          acc_rate < acc_min && 
          no_improve >= patience
      ) {
        std::cout << 
          "Stopping early at iteration " 
          << (z + 1)
          << " due to combined criterion: " 
          << "low acceptance AND no improvement.\n";
              
        early_stop_iter = (z + 1);
        stop_now = true;
    }

    if (lam_enabled && !stop_now) {
      
      const double hat_up =
        (uphill_attempted_in_win > 0)
      ? static_cast<double>(
          uphill_accepted_in_win
      ) / uphill_attempted_in_win
      : (
          (attempted_in_win > 0)
      ? static_cast<double>(
          accepted_in_win
      ) / attempted_in_win
      : 0.0
      );
        
      const double progress =
        static_cast<double>(z + 1) /
          static_cast<double>(n_iterations);
        
      const double a_star =
        target_acceptance(progress);
        
      const double num =
        std::log(clamp01(hat_up));
        
      const double den =
        std::log(clamp01(a_star));
        
      double ratio =
        (den != 0.0)
        ? (num / den)
        : 1.0;
        
      if (
          !(ratio > 0.0) ||
            !std::isfinite(ratio)
      ) {
        ratio = 1.0;
      }
        
      temp_curr *=
        std::pow(
          ratio,
          1.0 / lam_p
        );
        
      if (
          !std::isfinite(temp_curr) ||
            temp_curr < 1e-9
      ) {
        temp_curr = 1e-9;
      }
        
      if (temp_curr > 1e16) {
        temp_curr = 1e16;
      }
    }
      
    // Reset window counters
    attempted_in_win = 0;
    accepted_in_win  = 0;
      
    uphill_attempted_in_win = 0;
    uphill_accepted_in_win  = 0;
    }

    // -------- Write intermediate results to file --------
    if (
        trace_every > 0 &&
          trace_out.is_open() &&
          ((z + 1) % trace_every == 0)
    ) {
      
      const double overall_acc =
        (attempted_total > 0)
      ? static_cast<double>(
          accepted_total
      ) / attempted_total
      : std::numeric_limits<double>::quiet_NaN();
      
      const double T_curr =
        lam_enabled
        ? temp_curr
      : (
          temp /
            (
                1.0 +
                  cooling_rate_c * z
            )
      );
      
      const double avg_ms =
        (iter_count > 0)
        ? static_cast<double>(
            iter_ms_total
        ) / iter_count
      : 0.0;
      
      trace_out
      << (z + 1) << ","
      << best_score << ","
      << curr_eval << ","
      << T_curr << ","
      << accepted_total << ","
      << attempted_total << ","
      << overall_acc << ","
      << avg_ms
      << '\n';
      
      trace_out.flush();
    }
    
    // Checkpoint
    if (
        checkpoint_every > 0 &&
          !checkpoint_file.empty() &&
          (
              ((z + 1) % checkpoint_every == 0) ||
                stop_now
          )
    ) {
      
      SACheckpoint cp;
      
      // Metadata
      cp.dim_x = dim_x;
      cp.dim_y = dim_y;
      cp.n_actions = n_actions;
      cp.n_species = n_species;
      cp.n_habitats = n_habitats;
      cp.n_iterations = n_iterations;
      
      // Objective scaling
      cp.scale_cost   = scale_cost;
      cp.scale_config = scale_config;
      cp.scale_target = scale_target;
      
      // Main state
      cp.iteration = z + 1;
      cp.curr          = curr;
      cp.curr_eval     = curr_eval;
      cp.improve_count = improve_count;
      cp.best       = best;
      cp.best_score = best_score;
      cp.best_Fx    = best_Fx;
      cp.g_best         = g_best;
      cp.g_create_best  = g_create_best;
      cp.g_improve_best = g_improve_best;
      
      // Counters
      cp.attempted_total = attempted_total;
      cp.accepted_total  = accepted_total;
      cp.attempted_in_win = attempted_in_win;
      cp.accepted_in_win  = accepted_in_win;
      cp.no_improve = no_improve;
      cp.last_best  = last_best;
      cp.uphill_attempted_in_win = uphill_attempted_in_win;
      cp.uphill_accepted_in_win = uphill_accepted_in_win;
      cp.temp_curr = temp_curr;
      
      // Timings
      cp.iter_ms_total = iter_ms_total;
      cp.iter_count    = iter_count;
      
      // Histories
      cp.H_history  = H_history;
      cp.F_history  = F_history;
      cp.F1_history = F1_history;
      cp.F2_history = F2_history;
      cp.acc_rate_trace = acc_rate_trace;

      // RNG
      std::ostringstream rng_out;
      
      rng_out
      << rng
      << '\n'
      << runif;
      
      cp.rng_state =
      rng_out.str();
      
      write_sa_checkpoint(
        checkpoint_file,
        cp
      );
    }
    
    // Verbose output
    if (
        verbose &&
          ((z + 1) % win == 0)
    ) {
      
      std::cout
      << "[SA] it=" << (z + 1)
      << " T="
      << (
      lam_enabled
      ? temp_curr
      : (
          temp /
            (
                1.0 +
                  cooling_rate_c * z
            )
      )
      )
      << " Best H(x) = "
      << best_score
      << " Best F(x) = "
      << best_Fx
      << '\n'
      << std::flush;
    }
    
    
    // Only stop after state has been written
    if (stop_now) {
      break;
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
  out.diag.iterations_completed  = iterations_completed;
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
