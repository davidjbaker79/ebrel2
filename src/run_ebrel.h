//-------------------------- Run Ebrel model -----------------------------------

#ifndef RUN_EBREL_H
#define RUN_EBREL_H

#include "calibrate_objective_scaling.h"
#include <vector>
#include <stdexcept>
#include <limits>
#include <cstdint>
#include <string>

#include "offsets.h"
#include "species_plan.h"

//---------------------------- Main function -----------------------------------

// User supplied options for running ebrel model
struct RunEBRELOptions {
  double sigma            = 0.05;
  const std::vector<int8_t>* X0 = nullptr;   // X0 supply optional
  double base_prob_X0     = 0.95;
  double step_proportion  = 0.01;
  double step_probability = 0.05;
  int n_iterations        = 10000;
  double temp             = 2000;
  double cooling_rate_c   = 1;
  bool   lam_enabled      = false; // turn Lam-style online control on/off
  double lam_target_mid   = 0.44;  // target uphill acceptance during early/mid run
  double lam_target_final = 0.05;  // target uphill acceptance at the end
  double lam_hold_frac    = 0.60;  // fraction of run to hold lam_target_mid before decaying
  double lam_p            = 2.0;   // damping exponent in Ben-Ameur correction (>=1)
  int min_iterations      = 1000;  // require at least this many iterations
  int acceptance_window   = 1000;  // window length for acceptance rate
  double acceptance_thres = 0.01;  // "low" acceptance threshold
  int iter_no_improve     = 1000;  // consecutive iters with no meaningful improvement
  double improve_eps      = 1e-6;  // relative improvement needed to reset patience
  int rng_seed            = -1;
  int  write_every        = 0;     // 0 = disabled
  std::string trace_file  = ""; // empty = disabled
  bool verbose            = false;
};

// Input data structure
struct RunEBRELInput {
  int dim_x = 0;
  int dim_y = 0;
  int n_actions = 0;
  int n_species = 0;
  int n_habitats = 0;
  int universal_disp_thres = 20;
  int max_disp_steps = 10;
  int roi_cap  = 100;
  int cluster_gap_cells = 25;
  double alpha = 1.0;
  double beta  = 1.0;
  double gamma = 1.0;
  ObjectiveScales objective_scales;
  std::vector<int8_t> E;
  std::vector<uint8_t> U;
  std::vector<double> C, O, SD, W, SxH;
  std::vector<int>    D;
  std::vector<int8_t> X0;
  std::vector<uint8_t> LM;
  std::vector<int16_t> row_first_land;
  std::vector<int16_t> row_last_land;
  std::vector<int16_t> col_first_land;
  std::vector<int16_t> col_last_land;
  std::vector<int> cell_r, cell_c;
  std::vector<std::vector<std::size_t>> Etiles_per_h;
  OffsetsCache offsets_cache;
  RowRunsCache rowruns_cache;
  SpeciesPlan species_plan;
  std::vector<int> improve_count;
  std::vector<float> action_weight;
  double improve_w_default;
};

// Results structure
struct RunEBRELResult {
  std::vector<int8_t> X_best;
  double H_best = 0.0;
  std::vector<double> H_trace;
  std::vector<double> F_trace;
  std::vector<double> F1_trace;
  std::vector<double> F2_trace;
  std::vector<double> g_best;
  std::vector<double> g_create_best;
  std::vector<double> g_improve_best;
  int iterations_run = 0;
  // --- Diagnostics from SA ---
  std::vector<double> acc_rate_trace;  // acceptance per window
  int    early_stop_iter = -1;         // iteration where combined stop fired, or -1
  int    proposals = 0;                // total evaluated proposals (attempted_total)
  int    accepted  = 0;                // total accepted proposals
  double overall_acc = std::numeric_limits<double>::quiet_NaN(); // accepted/proposals
  // --- Timings ---
  long long iter_ms_total = 0;
  int iter_count = 0;
};

// Declare the wrapper
RunEBRELResult run_ebrel(const RunEBRELInput& in, const RunEBRELOptions& opt);

#endif // RUN_EBREL_H
