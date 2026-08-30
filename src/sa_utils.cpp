//---------------------Simulated annealing utility functions -------------------

#include "sa_utils.h"
#include "objective_utils.h"
#include "optimisation_utils.h"
#include "update_candidate.h"
#include "dispersal_utils.h"


#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <cstddef>

//---------------------- Local helper functions --------------------------------

namespace {

// Clamp funciton
inline double clamp01(double x, double eps = 1e-15) {
  if (x < eps) return eps;
  if (x > 1.0 - eps) return 1.0 - eps;
  return x;
}
  
  // For chi hat calculations
  inline double log_sum_exp_negE_over_T(const std::vector<double>& E, double invT) {
    if (E.empty()) return -std::numeric_limits<double>::infinity();
    
    double m = -std::numeric_limits<double>::infinity();
    for (double e : E) {
      const double v = -e * invT;
      if (v > m) m = v;
    }
    
    double s = 0.0;
    for (double e : E) {
      s += std::exp((-e * invT) - m);
    }
    
    if (s <= 0.0 || !std::isfinite(s)) {
      return -std::numeric_limits<double>::infinity();
    }
    
    return m + std::log(s);
  }
  
}

//--------------------- Main functions -----------------------------------------

// ---- Function to estimate initial temperature ----
InitTempResult estimate_initial_temperature_benameur_cpp(
    std::vector<int8_t> X_seed,
    const std::vector<double>& W,
    const std::vector<uint8_t>& U,
    const std::vector<double>& C,
    const std::vector<int8_t>& E,
    const std::vector<std::vector<std::size_t>>& Etiles_per_h,
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c,
    const RowRunsCache& rowruns_cache,
    const SpeciesPlan& species_plan,
    int   universal_disp_thres,
    int   max_disp_steps,
    int   roi_cap,
    const std::vector<uint8_t>& LM,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
    int   n_actions,
    int   n_species,
    int   n_habitats,
    int   dim_x,
    int   dim_y,
    double alpha,
    double beta,
    double gamma,
    double scale_cost,
    double scale_config,
    double scale_target,
    double step_proportion,
    double step_probability,
    double improve_total_weight, // default = 1
    int    num_samples,
    double chi0,                // target acceptance
    double p,                   // damping exponent in Ben-Ameur’s update (≥1)
    double tol_logchi,          // tolerance on |log χ̂ − log χ0|
    int    max_iters,           // cap on fixed-point iterations
    double T1,                  // optional initial guess; ≤0 => use heuristic
    int    max_tries_factor,    // limit on total neighbour proposals during
    int    rng_seed,
    bool   verbose
) {
  
  // ---- Checks ----
  if (!(chi0 > 0.0 && chi0 < 1.0)) throw std::invalid_argument("chi0 must be in (0,1)");
  if (!(p >= 1.0)) throw std::invalid_argument("p must be >= 1");
  if (num_samples <= 0) throw std::invalid_argument("num_samples must be > 0");
  if (dim_x <= 0 || dim_y <= 0) throw std::invalid_argument("dim_x and dim_y must be positive");
  if (n_actions <= 0 || n_species <= 0) throw std::invalid_argument("n_actions and n_species must be positive");
  if (!(step_proportion >= 0.0 && step_proportion <= 1.0)) throw std::invalid_argument("step_proportion must be in [0,1]");
  if (!(step_probability >= 0.0 && step_probability <= 1.0)) throw std::invalid_argument("step_probability must be in [0,1]");
  if (max_tries_factor <= 0) throw std::invalid_argument("max_tries_factor must be > 0");
  if (max_iters < 0) throw std::invalid_argument("max_iters must be >= 0");
  
  const std::size_t n_cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  if (n_cells > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("dim_x*dim_y too large for int indexing");
  }
  
  if (X_seed.size() != n_cells) {
    throw std::invalid_argument("X_seed must have length dim_x*dim_y");
  }
  
  const std::size_t Wsz = n_cells * static_cast<std::size_t>(n_habitats);
  const std::size_t Usz = n_cells * static_cast<std::size_t>(n_actions);
  if (W.size() != Wsz) throw std::invalid_argument("W must have length n_habitats*dim_x*dim_y");
  if (U.size() != Usz) throw std::invalid_argument("U must have length n_actions*dim_x*dim_y");
  
  // ---- RNG ----
  std::mt19937 rng;
  if (rng_seed >= 0) rng.seed(static_cast<uint32_t>(rng_seed));
  else { std::random_device rd; rng.seed(rd()); }
  
  std::uniform_real_distribution<> runif(0.0, 1.0);
  
  // ---- Vectors for storing uphill transitions ----
  std::vector<double> Emins, Emaxs, dpos;
  Emins.reserve(num_samples);
  Emaxs.reserve(num_samples);
  dpos.reserve(num_samples);
  
  // Improvement count
  std::vector<int> improve_count = initialise_improve_count(
    X_seed, 
    species_plan
  );
  
  // Initial
  double E_base = compute_H(
    X_seed, C, 
    alpha,
    beta,
    gamma,
    scale_cost,
    scale_config,
    scale_target,
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
  ).H;
  
  const int max_tries = std::max(num_samples * max_tries_factor, 100);
  int tries = 0;
  
  // Tune initial temp
  while (static_cast<int>(Emins.size()) < num_samples && tries++ < max_tries) {
    const int step_seed = static_cast<int>(rng());
    CandidateUpdate update = update_candidate(
      W, U, X_seed,
      step_proportion,
      step_probability,
      improve_total_weight,
      n_actions, 
      n_habitats,
      dim_x, dim_y,
      step_seed
    );
    
    if (update.changes.empty()) continue;
    const std::vector<int8_t>& candidate = update.X;
    
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
    
    // Compute H
    double E_cand = compute_H(
      candidate, C,
      alpha,
      beta,
      gamma,
      scale_cost,
      scale_config,
      scale_target,
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
    ).H;
    
    const double dE = E_cand - E_base;
    if (dE > 0.0) {
      Emins.push_back(E_base);
      Emaxs.push_back(E_cand);
      dpos.push_back(dE);
    }
    
    // Random walk — restore improve_count if not accepting
    if (runif(rng) < 0.5) {
      X_seed = candidate;
      E_base = E_cand;
      improve_count = std::move(improve_count_candidate);
    }
  }
  
  InitTempResult res{};
  res.diag.tries = tries;
  
  if (Emins.empty()) {
    // No uphill moves observed -> tiny T0 and diagnostics
    res.T0                 = 1e-6;
    res.diag.samples_used  = 0;
    res.diag.chi0          = chi0;
    res.diag.chi_hat_final = 0.0;
    res.diag.iters         = 0;
    res.diag.dE_mean       = std::numeric_limits<double>::quiet_NaN();
    res.diag.dE_median     = std::numeric_limits<double>::quiet_NaN();
    if (verbose) std::cerr << "[init_temp] No uphill transitions found; returning T0=1e-6\n";
    return res;
  }
  
  // stats
  const double dE_mean = std::accumulate(dpos.begin(), dpos.end(), 0.0) / (double)dpos.size();
  std::nth_element(dpos.begin(), dpos.begin() + dpos.size()/2, dpos.end());
  const double dE_median = dpos[dpos.size()/2];
  
  // Initial temp guess
  double T = T1;
  if (!(T > 0.0)) {
    T = dE_median / (-std::log(chi0));
    if (!std::isfinite(T) || T <= 0.0) T = 1.0;
  }
  
  // chi_hat(T)
  auto chi_hat = [&](double Tcur) {
    if (!(Tcur > 0.0) || !std::isfinite(Tcur)) return clamp01(0.0);
    const double invT = 1.0 / Tcur;
    const double lnum = log_sum_exp_negE_over_T(Emaxs, invT);
    const double lden = log_sum_exp_negE_over_T(Emins, invT);
    const double val  = std::exp(lnum - lden);
    return clamp01(val);
  };
  
  // Ben-Ameur iterations
  const double log_chi0 = std::log(clamp01(chi0));
  int iters = 0;
  for (; iters < max_iters; ++iters) {
    const double ch  = chi_hat(T);
    const double lch = std::log(clamp01(ch));
    
    if (std::fabs(lch - log_chi0) <= tol_logchi) break;
    
    double ratio = (log_chi0 != 0.0) ? (lch / log_chi0) : 1.0; // both negative
    if (!(ratio > 0.0) || !std::isfinite(ratio)) ratio = 1.0;
    T *= std::pow(ratio, 1.0 / p);
    
    if (!std::isfinite(T) || T <= 1e-16) T = 1e-6;
    if (T > 1e16) T = 1e16;
  }
  
  res.T0 = T;
  res.diag.iters         = iters;
  res.diag.samples_used  = static_cast<int>(Emins.size());
  res.diag.chi0          = chi0;
  res.diag.chi_hat_final = chi_hat(T);
  res.diag.dE_mean       = dE_mean;
  res.diag.dE_median     = dE_median;
  
  if (verbose) {
    std::cout << "[init_temp] T0=" << res.T0
              << " | chi_hat(T0)=" << res.diag.chi_hat_final
              << " | samples=" << res.diag.samples_used
              << " | iters=" << res.diag.iters
              << " | dE_med=" << res.diag.dE_median
              << " | dE_mean=" << res.diag.dE_mean
              << std::endl;
  }
  return res;
}



