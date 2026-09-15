// -------------------------------Simulated annealing----------------------------

#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <string>

#include "offsets.h"
#include "species_plan.h"

// ------------------------------- Main functions ------------------------------

// For SA diagnostics
struct SADiagnostics {
  std::vector<double> acceptance_rate_trace; // acceptance per 'acceptance_window'
  int                 early_stop_iter = -1;  // iteration (1-based) is stopped early, -1 otherwise
  int                 attempted_total = 0;   // total proposals attempted (non no-op)
  int                 accepted_total  = 0;   // total accepted
  int                 iterations_completed = 0; // number of iterations
  long long           iter_ms_total   = 0;   // total iterations time in ms
  int                 iter_count      = 0;   // number of iterations
};

// Result container for simulated annealing (Rcpp-free)
struct SAResult {
  std::vector<int8_t> X_best;   // flattened [n_cells]
  double H_best;                // best objective
  std::vector<double> H_trace;  // candidate H evaluations
  std::vector<double> F_trace;  // optional Fx trace from compute_H
  std::vector<double> F1_trace; // optional F1 trace
  std::vector<double> F2_trace; // optional F2 trace
  std::vector<double> g_best;   // Short-fall in targets
  std::vector<double> g_create_best;   // Creation gain
  std::vector<double> g_improve_best;   // Improvement gain
  SADiagnostics diag;
};

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
    // ---------- User specified 
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
);

struct SACheckpoint {
  
  uint32_t version = 1;
  
  int dim_x = 0;
  int dim_y = 0;
  int n_actions = 0;
  int n_species = 0;
  int n_habitats = 0;
  int n_iterations = 0;
  
  // Number of completed SA iterations
  int iteration = 0;
  
  // Current chain state
  std::vector<int8_t> curr;
  double curr_eval = 0.0;
  std::vector<int> improve_count;
  
  // Best encountered state
  std::vector<int8_t> best;
  double best_score = 0.0;
  double best_Fx = 0.0;
  
  std::vector<double> g_best;
  std::vector<double> g_create_best;
  std::vector<double> g_improve_best;
  
  // Acceptance / stopping state
  int attempted_total = 0;
  int accepted_total = 0;
  
  int attempted_in_win = 0;
  int accepted_in_win = 0;
  
  int no_improve = 0;
  double last_best = 0.0;
  
  // Lam state
  int uphill_attempted_in_win = 0;
  int uphill_accepted_in_win = 0;
  double temp_curr = 0.0;
  
  // Timings
  long long iter_ms_total = 0;
  int iter_count = 0;
  
  // Histories
  std::vector<double> H_history;
  std::vector<double> F_history;
  std::vector<double> F1_history;
  std::vector<double> F2_history;
  std::vector<double> acc_rate_trace;
  
  // Serialized RNG state
  std::string rng_state;
};


void write_sa_checkpoint(
    const std::string& checkpoint_file,
    const SACheckpoint& cp
);

SACheckpoint read_sa_checkpoint(
    const std::string& checkpoint_file
);

#endif // SIMULATED_ANNEALING_H