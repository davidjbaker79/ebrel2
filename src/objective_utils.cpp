
//-------------------------- Ebrel objective functions -------------------------

#include "objective_utils.h"
#include "dispersal_utils.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>  // std::runtime_error, std::invalid_argument
#include <string>     // std::to_strings
//#include <iostream>


//-------------------------- Local helper functions ----------------------------
namespace {

  // Get number of cells
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  inline std::size_t suitable_idx(const SpeciesPlan& P, int sp, int h) {
    return static_cast<std::size_t>(sp) * static_cast<std::size_t>(P.n_habitats) + static_cast<std::size_t>(h);
  }

  inline bool is_suitable(const SpeciesPlan& P, int sp, int h) {
    return P.suitable[suitable_idx(P, sp, h)] != 0u;
  }

  // For Euclidean distance
  inline int row_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile / static_cast<std::size_t>(dim_x));
  }

  // For Euclidean distance
  inline int col_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile % static_cast<std::size_t>(dim_x));
  }

  // Compute discounted seeds
  inline double compute_m_i_from_plan(
      int sp,
      const SpeciesPlan& species_plan,
      const std::vector<int8_t>& X,
      const std::vector<int8_t>& E)
  {
    double m_i = 0.0;
    
    const double omega =
      species_plan.marginal_occupancy_weight;
    
    const uint32_t start =
      species_plan.seed_start[static_cast<std::size_t>(sp)];
    
    const uint32_t len =
      species_plan.seed_len[static_cast<std::size_t>(sp)];
    
    for (uint32_t j = 0; j < len; ++j) {
      
      const std::size_t seed_idx =
        static_cast<std::size_t>(start + j);
      
      const std::size_t k =
        static_cast<std::size_t>(
          species_plan.seed_pool[seed_idx]
        );
      
      const double baseline_w =
        static_cast<double>(
          species_plan.seed_weight[seed_idx]
        );
      
      const int he = E[k];
      const int hx = X[k];
      
      const bool existing_suitable =
        he >= 0 &&
        he < species_plan.n_habitats &&
        is_suitable(species_plan, sp, he);
      
      if (hx < 0) {
        
        // No action: retain baseline occupancy.
        m_i += baseline_w;
        
      } else if (hx < species_plan.n_habitats) {
        
        const bool new_suitable =
          is_suitable(species_plan, sp, hx);
        
        if (existing_suitable) {
          
          if (new_suitable) {
            // Suitable main habitat replaced by another suitable habitat:
            // cell remains fully occupied.
            m_i += 1.0;
          } else {
            // Suitable main habitat lost, but marginal/sub-grid habitat remains.
            m_i += omega;
          }
          
        } else {
          
          // Species was initially supported only by marginal/sub-grid habitat.
          // That marginal occupancy survives conversion of the mapped main habitat.
          m_i += omega;
        }
        
      } else {
        
        // Improvement: retain baseline occupancy.
        m_i += baseline_w;
      }
    }
    
    return m_i;
  }

  // Sum of Euclidean distances for unordered pairs within one set (i<j).
  // - can multiply by 2 to obtain the ordered-pair sum to match JAE
  double sum_pairwise_self_unordered(const std::vector<std::size_t>& tiles, int dim_x) {
    const std::size_t n = tiles.size();
    if (n < 2) return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      const int ri = row_of(tiles[i], dim_x), ci = col_of(tiles[i], dim_x);
      for (std::size_t j = i + 1; j < n; ++j) {
        const int rj = row_of(tiles[j], dim_x), cj = col_of(tiles[j], dim_x);
        const double dr = static_cast<double>(ri - rj);
        const double dc = static_cast<double>(ci - cj);
        s += std::sqrt(dr*dr + dc*dc);
      }
    }
    return s;
  }

  // Sum of Euclidean distances across two sets (e.g. E × X0).
  double sum_pairwise_cross(const std::vector<std::size_t>& A,
                            const std::vector<std::size_t>& B,
                            int dim_x) {
    if (A.empty() || B.empty()) return 0.0;
    double s = 0.0;
    for (std::size_t ia = 0; ia < A.size(); ++ia) {
      const int ra = row_of(A[ia], dim_x), ca = col_of(A[ia], dim_x);
      for (std::size_t ib = 0; ib < B.size(); ++ib) {
        const int rb = row_of(B[ib], dim_x), cb = col_of(B[ib], dim_x);
        const double dr = static_cast<double>(ra - rb);
        const double dc = static_cast<double>(ca - cb);
        s += std::sqrt(dr*dr + dc*dc);
      }
    }
    return s;
  }

}

// -------------------- Objective functions ------------------------------------

// F1: sum_{h,cell} X * C  (habitat-major layout: idx = h*cells + tile)
double compute_F1(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  int dim_x, 
                  int dim_y) {

  const int cells = n_cells(dim_x, dim_y);
  double f1 = 0.0;

  for (int cell = 0; cell < cells; ++cell) {
    const int h = X[cell];
    if (h >= 0) {
      const std::size_t idx =
        static_cast<std::size_t>(h) * static_cast<std::size_t>(cells)
      + static_cast<std::size_t>(cell);
      f1 += C[idx];
    }
  }

  return f1;
}

// F2: total = w_xx * (X–X sum once) + w_ex * (E–X sum once)
//  - original formulation used ordered pairwise sum, but this weights xx more than xe
//  - in this version, where creation and improvement are allowed
double compute_F2(const std::vector<int8_t>& X,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  int n_habitats, int dim_x, int dim_y) {

  // Number of cells
  const int cells = n_cells(dim_x, dim_y);
  
  // Build X tiles per habitat in one pass over cells
  std::vector<std::vector<std::size_t>> Xtiles_per_h(static_cast<std::size_t>(n_habitats));
  // Optional: rough reserve heuristic (tune or drop)
  for (int h = 0; h < n_habitats; ++h) Xtiles_per_h[static_cast<std::size_t>(h)].reserve(cells / 16);

  for (int tile = 0; tile < cells; ++tile) {
    const int h = X[static_cast<std::size_t>(tile)];
    if (h >= 0 && h < n_habitats) { // exclude improvement action (h == n_habitats)
      Xtiles_per_h[static_cast<std::size_t>(h)].push_back(
          static_cast<std::size_t>(tile));
    }
  }

  double total = 0.0;

  for (int h = 0; h < n_habitats; ++h) {
    const auto& Etiles = Etiles_per_h[static_cast<std::size_t>(h)];
    const auto& Xtiles = Xtiles_per_h[static_cast<std::size_t>(h)];

    if (!Xtiles.empty()) {
      total += sum_pairwise_self_unordered(Xtiles, dim_x);     // X–X
      if (!Etiles.empty()) total += sum_pairwise_cross(Etiles, Xtiles, dim_x); // E–X
    }
  }

  return total;
}

// ---- Compute H
HResult compute_H(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  double alpha,
                  double beta,
                  double gamma,
                  double delta,
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
                  const std::vector<int>& improve_count  // n_s: improved cells per species (changes each SA step)
){

  // ---- F1 and F2 -----
  const double f1 = compute_F1(X, C, dim_x, dim_y);
  const double f2 = compute_F2(X, Etiles_per_h, n_habitats, dim_x, dim_y);

  // ---- Compute G ----
  std::vector<double> G = compute_G(
    X,
    LM,
    row_first_land, row_last_land,
    col_first_land, col_last_land,
    cell_r, cell_c,
    E,
    dim_x, dim_y,
    universal_disp_thres, max_disp_steps, roi_cap,
    rowruns_cache,
    species_plan);

  // ---- Shortfall per species relative to target o_i * m_i with creation and improvement ----
  std::vector<double> g(static_cast<std::size_t>(n_species), 0.0);
  std::vector<double> g_improve(static_cast<std::size_t>(n_species), 0.0);
  std::vector<double> g_create(static_cast<std::size_t>(n_species), 0.0);
  
  // ---- Missed targets
  int n_target_species = 0;
  int n_targets_missed = 0;
  
  constexpr double target_tol = 1e-8;
  
  for (int sp = 0; sp < n_species; ++sp) {

    const std::size_t spi = static_cast<std::size_t>(sp);
    
    if (!species_plan.active[spi]) {
      g[spi] = 0.0;
      g_create[spi] = 0.0;
      g_improve[spi] = 0.0;
      continue;
    }
    
    const double target = species_plan.O_n[spi];
    
    if (target <= 0.0) {
      g[spi] = 0.0;
      g_create[spi] = 0.0;
      g_improve[spi] = 0.0;
      continue;
    }
    
    // Species has a genuine target in this optimisation.
    ++n_target_species;
    
    const double m0 =
    species_plan.seed_weight_sum[spi];
    
    const double kept =
      compute_m_i_from_plan(
        sp,
        species_plan,
        X,
        E
      );
    
    const double lost = m0 - kept;
    
    // G is already habitat-weighted inside compute_G().
    // Express creation as a proportion of this species' target.
    g_create[spi] = (G[spi] - lost) / target;
    
    // Improvement contributes in the same absolute species-target units as
    // creation. action_weight for the improvement action is the contribution
    // of one selected, suitable and reachable improvement cell.
    const int improve_action = n_actions - 1;
    
    const double improve_benefit =
      static_cast<double>(improve_count[spi]) *
      static_cast<double>(
        species_plan.action_weight[
          static_cast<std::size_t>(improve_action)
        ]
      );
    
    // Express that absolute benefit as a proportion of this species' target.
    g_improve[spi] = improve_benefit / target;
    
    // Residual target shortfall: 1 = no progress; 0 = target met.
    const double achieved = g_create[spi] + g_improve[spi];
    
    g[spi] = std::max(0.0, 1.0 - achieved);

    // Count species that have not achieved their target.
    if (g[spi] > target_tol) {
      ++n_targets_missed;
    }
    
  }

  const double gx_val = std::accumulate(g.begin(), g.end(), 0.0);
  
  const double missed_target_prop =
    (n_target_species > 0)
    ? static_cast<double>(n_targets_missed) /
      static_cast<double>(n_target_species)
    : 0.0;
  
  // The three raw components are standardised using fixed calibration
  // values estimated before temperature tuning and simulated annealing.
  //
  // User preferences are then applied to those standardised components.
  // For example, beta = 1 and gamma = 10 means that a one-reference-unit
  // reduction in target shortfall is weighted ten times more strongly than
  // a one-reference-unit reduction in configuration penalty.
  if (!std::isfinite(scale_cost) || scale_cost <= 0.0) {
    throw std::invalid_argument(
        "compute_H: scale_cost must be finite and positive"
    );
  }
  
  if (!std::isfinite(scale_config) || scale_config <= 0.0) {
    throw std::invalid_argument(
        "compute_H: scale_config must be finite and positive"
    );
  }
  
  if (!std::isfinite(scale_target) || scale_target <= 0.0) {
    throw std::invalid_argument(
        "compute_H: scale_target must be finite and positive"
    );
  }
  
  const double F1_standardised = f1 / scale_cost;
  const double F2_standardised = f2 / scale_config;
  const double g_standardised  = gx_val / scale_target;
  
  const double Fx_val =
    alpha * F1_standardised +
    beta  * F2_standardised;
  
  const double H_val =
    Fx_val +
    gamma * g_standardised +
    delta * missed_target_prop;

  return { 
      H_val, 
      Fx_val, 
      gx_val, 
      f1, 
      f2, 
      G, 
      g, 
      g_create, 
      g_improve,
      missed_target_prop,
      n_targets_missed,
      n_target_species
    };
}
