//------------------------ Objective function utilities ------------------------

#ifndef OBJECTIVE_UTILS_H
#define OBJECTIVE_UTILS_H

#include <vector>
#include <cmath>
#include "dispersal_utils.h"  // compute_G

//-------------------------- Main Functions ------------------------------------

// Return structure for H(x)
struct HResult {
  double H;                       // Full weighted, standardised objective
  double Fx;                      // Weighted standardised F1 + F2 contribution
  double gx;                      // Raw summed residual target shortfall
  double F1;                      // Raw cost component
  double F2;                      // Raw configuration component
  std::vector<double> G;          // Raw creation benefit by species
  std::vector<double> g;          // Residual target shortfall by species
  std::vector<double> g_create;   // Creation contribution by species
  std::vector<double> g_improve;  // Improvement contribution by species
};

// Declare F1 objective function
double compute_F1(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  int dim_x,
                  int dim_y);

// Declare F2 objective function
double compute_F2(const std::vector<int8_t>& X,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  int n_habitats,
                  int dim_x,
                  int dim_y);

// Euclidean distance helper functions
double euclidean_distance(const std::vector<double>& a,
                          const std::vector<double>& b);

// Pairwise distances helper functions
double pairwise_distances(const std::vector<std::vector<double>>& A,
                          const std::vector<std::vector<double>>& B);

// Declare function to compute H(x)
HResult compute_H(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  double alpha,
                  double beta,
                  double gamma,
                  double scale_cost,
                  double scale_config,
                  double scale_target,
                  int n_actions,
                  int n_species,
                  int n_habitats,
                  int dim_x,
                  int dim_y,
                  int universal_disp_thres,
                  int max_disp_steps,
                  int roi_cap,
                  const std::vector<uint8_t>& LM,
                  const std::vector<int16_t>& row_first_land,
                  const std::vector<int16_t>& row_last_land,
                  const std::vector<int16_t>& col_first_land,
                  const std::vector<int16_t>& col_last_land,
                  const std::vector<int8_t>& E,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  const std::vector<int>& cell_r,
                  const std::vector<int>& cell_c,
                  const RowRunsCache& rowruns_cache,
                  const SpeciesPlan& species_plan,
                  std::vector<int> improve_count
                  );



#endif // OBJECTIVE_UTILS_H
