//-------------------------- Simulated Annealing Utilities ---------------------

#ifndef SA_UTILS_H
#define SA_UTILS_H

#include <vector>
#include <cstdint>
#include "dispersal_utils.h"

// Diagnostics returned with the temperature estimate
struct InitTempDiag {
  double chi_hat_final;   // final estimated acceptance at T0
  double chi0;            // target acceptance
  int    iters;           // Ben-Ameur iterations used
  int    samples_used;    // number of uphill transitions collected
  double dE_median;       // median positive deltaE in samples
  double dE_mean;         // mean positive deltaE in samples
  int    tries;           // attempts made to gather samples
};

// Result bundle
struct InitTempResult {
  double T0;
  InitTempDiag diag;
};

//--------------------------------- Main functions -----------------------------


/**
 * Estimate initial SA temperature T0 targeting acceptance probability chi0,
 * using Ben-Ameur (2004).
 *
 * - Uses your compute_H(...) and update_candidate(...) to sample uphill moves.
 * - W, U, ... match your SA function signatures.
 *
 * step_proportion / step_probability are passed to update_candidate(...)
 * so the neighborhood used for sampling matches your SA settings.
 */
InitTempResult estimate_initial_temperature_benameur_cpp(
    std::vector<int8_t> X_seed,           // starting state (copy)
    const std::vector<double>& W,         // distance weights
    const std::vector<uint8_t>& U,
    const std::vector<double>& C,
    const std::vector<int8_t>& E,
    const std::vector<std::vector<std::size_t>>& Etiles_per_h,
    const std::vector<int>&    cell_r,
    const std::vector<int>&    cell_c,
    const RowRunsCache& rowruns_cache,
    const SpeciesPlan& species_plan,
    int    universal_disp_thres,
    int    max_disp_steps,
    int    roi_cap,
    const  std::vector<uint8_t>& LM,
    const  std::vector<int16_t>& row_first_land,
    const  std::vector<int16_t>& row_last_land,
    const  std::vector<int16_t>& col_first_land,
    const  std::vector<int16_t>& col_last_land,
    int    n_actions,
    int    n_species,
    int    m_habitats,
    int    dim_x,
    int    dim_y,
    double alpha,
    double beta,
    double gamma,
    double scale_cost,
    double scale_config,
    double scale_target,
    // neighbourhood controls (forwarded to update_candidate)
    double step_proportion  = 0.05,
    double step_probability = 0.15,
    // default value for per-cell selection prob for improvement   
    double improve_total_weight = 1,
    // estimator controls
    int    num_samples      = 400,
    double chi0             = 0.8,
    double p                = 2.0,
    double tol_logchi       = 1e-3,
    int    max_iters        = 50,
    double T1               = -1.0,   // <=0 => use heuristic
    int    max_tries_factor = 50,
    int    rng_seed         = -1,     // <0 => nondeterministic
    bool   verbose          = false
);

#endif // SA_UTILS_H
