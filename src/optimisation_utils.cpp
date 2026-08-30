#include "optimisation_utils.h"
#include <queue>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "species_plan.h"

//-------------------- Local helper functions ----------------------------------
namespace {

  // For storing coordinates
  struct Cell { int row, col; };

  // For checking within the model boundary
  inline bool in_bounds(int r, int c, int dim_y, int dim_x) {
    return r >= 0 && c >= 0 && r < dim_y && c < dim_x;
  }

}

//------------------- Main functions -------------------------------------------

// W contains creation-action proposal weights only.
// Creation action h corresponds to habitat h, for h = 0 ... n_actions - 1.
// The improvement action is n_habitats (= n_actions + 1) and is handled elsewhere.
std::vector<double> compute_distance_weights(
    const std::vector<int8_t>& E,
    const std::vector<uint8_t>& U,
    int n_habitats,
    int dim_x,
    int dim_y,
    double sigma
) {
  const int n_cells = dim_y * dim_x;
  std::vector<double> W(n_cells * n_habitats, 0.0);

  // 8-neighbourhood offsets
  const int dr[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
  const int dc[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

  auto idx = [n_cells](int cell, int h) {
    return h * n_cells + cell;
    };

  const int INF = std::numeric_limits<int>::max();

  for (int h = 0; h < n_habitats; ++h) {
    // Int hops for Chebyshev
    std::vector<int> dist(n_cells, INF);
    std::queue<Cell> q;

    // Seed with presence cells
    for (int cell = 0; cell < n_cells; ++cell) {
      if (E[cell] == h) {
        int r = cell / dim_x, c = cell % dim_x;
        dist[cell] = 0;
        q.push(Cell{r, c});
      }
    }

    // BFS with unit-cost edges (diag == ortho == 1)
    while (!q.empty()) {
      Cell curr = q.front(); q.pop();
      int curr_idx = curr.row * dim_x + curr.col;
      int curr_dist = dist[curr_idx];

      for (int k = 0; k < 8; ++k) {
        int nr = curr.row + dr[k], nc = curr.col + dc[k];
        if (!in_bounds(nr, nc, dim_y, dim_x)) continue;

        int nb = nr * dim_x + nc;
        int new_dist = curr_dist + 1;
        if (new_dist < dist[nb]) {
          dist[nb] = new_dist;
          q.push(Cell{nr, nc});
        }
      }
    }

    // Convert to weights, mask U and E
    double total_weight = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {

      bool unavailable = (U[idx(cell, h)] != 0) || (E[cell] == h);

      if (unavailable || dist[cell] == INF)  {
        W[idx(cell, h)] = 0.0;
      } else {
        double w = std::exp(-sigma * dist[cell]);
        W[idx(cell, h)] = w;
        total_weight += w;
      }
    }

    // Normalize column h to sum 1
    if (total_weight > 0.0) {
      for (int cell = 0; cell < n_cells; ++cell) {
        double& w = W[idx(cell, h)];
        if (w > 0.0) w /= total_weight;
      }
    }
  }

  return W;
}

// Initialise improvement count
std::vector<int> initialise_improve_count(
    const std::vector<int8_t>& X,
    const SpeciesPlan& P
) {
  const int n_species = P.n_species;
  const int improve_action = P.n_habitats;  // last action
  
  std::vector<int> improve_count(n_species, 0);
  
  for (int sp = 0; sp < n_species; ++sp) {
    const uint32_t start = P.seed_start[sp];
    const uint32_t len   = P.seed_len[sp];
    
    int count = 0;
    
    for (uint32_t i = 0; i < len; ++i) {
      const uint32_t k = P.seed_pool[start + i];
      
      if (X[k] == improve_action) {
        count++;
      }
    }
    
    improve_count[sp] = count;
  }
  
  return improve_count;
}

