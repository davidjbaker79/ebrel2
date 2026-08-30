//------------------------------ Dispersal-aware G function --------------------

#include "species_plan.h"
#include "dispersal_utils.h"
#include "offsets.h"

#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <limits> //  for std::numeric_limits
#include <cmath>   // for std::abs

#ifdef _OPENMP
  #include <omp.h>
#endif

//-------------------------- Local helpers -------------------------------------

namespace {

inline bool visited_test_and_set(std::size_t j,
                                 std::vector<uint64_t>& bits,
                                 std::vector<uint32_t>& word_stamp,
                                 uint32_t cur_stamp)
{
  const std::size_t w = j >> 6;              // /64
  const uint64_t mask = 1ULL << (j & 63);    // %64

  if (word_stamp[w] != cur_stamp) {
    word_stamp[w] = cur_stamp;
    bits[w] = 0ULL;
  }
  const uint64_t old = bits[w];
  if (old & mask) return true;
  bits[w] = old | mask;
  return false;
}


}

//------------------------- G Function -----------------------------------------

/* Compute G(x): per-species count of reachable newly created suitable habitat (X),
 * using E (existing) and suitable-X as transit
 * The index conventions are:
 * -1 = no action
 * 0..n_habitats-1 = creation habitat types
 */
std::vector<double> compute_G(
    const std::vector<int8_t>& X,
    const std::vector<uint8_t>& LM,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c,
    const std::vector<int8_t>& E,
    int dim_x, int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    const RowRunsCache& rowruns_cache,
    const SpeciesPlan& species_plan) 
  {

  const int n_species = species_plan.n_species;
  const int n_habitats = species_plan.n_habitats;

  std::vector<double> result(n_species, 0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  // Habitat, transit and hab weight ROIS-local
  std::vector<uint8_t> transit, habitat;
  std::vector<float> hab_weight;  // ROI-local, weight of habitat contribution
  
  // For tracking BFS
  std::vector<uint64_t> visited_bits;
  std::vector<uint32_t> visited_word_stamp;
  uint32_t cur_stamp = 1;
  std::vector<std::size_t> frontier, next;
  std::vector<std::vector<std::pair<int,int>>> row_intervals;
  std::vector<int> touched_rows;
  const std::size_t frontier_group_threshold = 512;

  // Species loop
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int sp = 0; sp < n_species; ++sp) {

    if (!species_plan.active[sp]) { result[sp] = 0.0; continue; }

    const int16_t disp = species_plan.disp[sp];
    const int r0 = species_plan.r0[sp];
    const int r1 = species_plan.r1[sp];
    const int c0 = species_plan.c0[sp];
    const int c1 = species_plan.c1[sp];
    const int W = species_plan.W[sp];
    const int H = species_plan.H[sp];
    const std::size_t Nroi = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);

    const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_habitats);
    const uint8_t* suitable_h_flag = species_plan.suitable.data() + base_h;

    // ---- Fast path (global dispersal) still uses ROI but no BFS
    if (species_plan.fast_path[sp]) {
      double total = 0.0;
      for (int r = r0; r <= r1; ++r) {
        std::size_t g =
          static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x) +
          static_cast<std::size_t>(c0);
        for (int c = c0; c <= c1; ++c, ++g) {
          if (!LM[g]) continue;
          const int hx = X[g];
          // Exclude improvement (n_habitats) and unassigned (-1) from G
          if (hx != -1 && hx != n_habitats && suitable_h_flag[hx]) 
            total += species_plan.action_weight[
          static_cast<std::size_t>(hx)];
        }
      }
      result[sp] = total;
      continue;
    }

    // Scratch ROI-local
    if (Nroi > transit.capacity()) { 
      transit.reserve(Nroi); 
      habitat.reserve(Nroi); 
      hab_weight.reserve(Nroi); 
    }
    transit.resize(Nroi);
    habitat.resize(Nroi);
    hab_weight.resize(Nroi);
    
    frontier.clear(); next.clear();
    frontier.reserve(std::min<std::size_t>(Nroi, 1u << 20));
    next.reserve(std::min<std::size_t>(Nroi, 1u << 20));

    const std::size_t nwords = (Nroi + 63) >> 6;
    if (nwords > visited_bits.capacity()) { 
      visited_bits.reserve(nwords); 
      visited_word_stamp.reserve(nwords); 
      }
    visited_bits.resize(nwords);
    if (visited_word_stamp.size() < nwords) visited_word_stamp.resize(nwords, 0u);

    ++cur_stamp;
    if (cur_stamp == 0u) { 
      std::fill(visited_word_stamp.begin(), visited_word_stamp.end(), 0u); 
      cur_stamp = 1u; 
      }

    if ((int)row_intervals.size() < H) row_intervals.resize(H);
    touched_rows.clear();
    touched_rows.reserve(std::min(H, 4096));

    // ---- Build transit/habitat/hab_weight masks
    std::size_t roiH = 0;
    for (int r = r0; r <= r1; ++r) {
      const std::size_t rowOff = static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W);

      if (row_last_land[r] < c0 || row_first_land[r] > c1) {
        std::fill_n(&habitat[rowOff], W, 0u);
        std::fill_n(&transit[rowOff], W, 0u);
        continue;
      }

      const int cc0 = std::max(c0, static_cast<int>(row_first_land[r]));
      const int cc1 = std::min(c1, static_cast<int>(row_last_land[r]));

      const std::size_t left = static_cast<std::size_t>(cc0 - c0);
      const std::size_t right_start = static_cast<std::size_t>(cc1 - c0 + 1);
      const std::size_t right_len = static_cast<std::size_t>(W) - right_start;

      std::fill_n(&habitat[rowOff], left, 0u);
      std::fill_n(&transit[rowOff], left, 0u);
      if (right_len > 0) {
        std::fill_n(&habitat[rowOff + right_start], right_len, 0u);
        std::fill_n(&transit[rowOff + right_start], right_len, 0u);
      }

      std::size_t g = static_cast<std::size_t>(r) * 
        static_cast<std::size_t>(dim_x) + static_cast<std::size_t>(cc0);
      std::size_t idx = rowOff + static_cast<std::size_t>(cc0 - c0);

      for (int c = cc0; c <= cc1; ++c, ++idx, ++g) {
        uint8_t hmask = 0u, tmask = 0u;
        float hw = 0.0f;
        
        const int hx = X[g];
        const bool is_improve = (hx == n_habitats); // improvement action is n_habitats in 0-indexed
        
        if (!is_improve && hx != -1 && suitable_h_flag[hx]) {
          // Creation actions (0..n_habitats-1): counts as new habitat AND transit
          hmask = 1u; 
          tmask = 1u;
          hw = static_cast<float>(species_plan.action_weight[static_cast<std::size_t>(hx)]);
        } else {
          // No action (-1) or improvement (n_habitats): transit from existing habitat only
          const int he = E[g];
          if (he != -1 && suitable_h_flag[he]) tmask = 1u;
          hw = 0.0f;
          // hmask stays 0
        }
        habitat[idx] = hmask;
        transit[idx] = tmask;
        hab_weight[idx] = hw;
        roiH += hmask;
      }
    }

    // No suitable habitat
    if (roiH == 0) { result[sp] = 0.0; continue; }

    // ---- Seed frontier from flat pool
    const uint32_t start = species_plan.seed_start[sp];
    const uint32_t len   = species_plan.seed_len[sp];

    for (uint32_t t = 0; t < len; ++t) {
      const std::size_t k = static_cast<std::size_t>(species_plan.seed_pool[start + t]);

      const int r = cell_r[k], c = cell_c[k];
      if (r < r0 || r > r1 || c < c0 || c > c1) continue;

      const std::size_t idx = static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W)
        + static_cast<std::size_t>(c - c0);

      if (!transit[idx]) continue;
      if (visited_test_and_set(idx, visited_bits, visited_word_stamp, cur_stamp)) continue;
      frontier.push_back(idx);
    }

    // Nowhere to go
    if (frontier.empty()) { result[sp] = 0.0; continue; }

    const RowRuns& runs = rowruns_cache[static_cast<std::size_t>(disp)];
    double g_sp = 0.0;

    for (int step = 0; !frontier.empty() && (max_disp_steps <= 0 || step < max_disp_steps); ++step) {
      next.clear();

      if (frontier.size() < frontier_group_threshold) {
        // small frontier path
        for (std::size_t i : frontier) {
          const int fr = static_cast<int>(i / static_cast<std::size_t>(W));
          const int fc = static_cast<int>(i - static_cast<std::size_t>(fr) * 
                                          static_cast<std::size_t>(W));

          for (const RowRun& rr : runs) {
            const int nr = fr + static_cast<int>(rr.dr);
            if (static_cast<unsigned>(nr) >= static_cast<unsigned>(H)) continue;

            int a = fc + static_cast<int>(rr.dc_min);
            int b = fc + static_cast<int>(rr.dc_max);
            if (b < 0 || a >= W) continue;
            if (a < 0) a = 0;
            if (b >= W) b = W - 1;

            const std::size_t base = static_cast<std::size_t>(nr) * static_cast<std::size_t>(W);
            for (int nc = a; nc <= b; ++nc) {
              const std::size_t j = base + static_cast<std::size_t>(nc);
              if (!transit[j]) continue;
              if (visited_test_and_set(j, visited_bits, visited_word_stamp, cur_stamp)) continue;
              if (habitat[j]) g_sp += static_cast<double>(hab_weight[j]);
              next.push_back(j);
            }
          }
        }
      } else {
        // grouped frontier path
        for (int rr : touched_rows) row_intervals[rr].clear();
        touched_rows.clear();

        for (std::size_t i : frontier) {
          const int fr = static_cast<int>(i / static_cast<std::size_t>(W));
          const int fc = static_cast<int>(i - static_cast<std::size_t>(fr) * 
                                          static_cast<std::size_t>(W));

          for (const RowRun& rr : runs) {
            const int nr = fr + static_cast<int>(rr.dr);
            if (static_cast<unsigned>(nr) >= static_cast<unsigned>(H)) continue;

            int a = fc + static_cast<int>(rr.dc_min);
            int b = fc + static_cast<int>(rr.dc_max);
            if (b < 0 || a >= W) continue;
            if (a < 0) a = 0;
            if (b >= W) b = W - 1;

            auto& vec = row_intervals[nr];
            if (vec.empty()) touched_rows.push_back(nr);
            vec.emplace_back(a, b);
          }
        }

        for (int nr : touched_rows) {
          auto& ints = row_intervals[nr];
          if (ints.empty()) continue;
          std::sort(ints.begin(), ints.end());

          int cur_a = ints[0].first;
          int cur_b = ints[0].second;

          const std::size_t base = static_cast<std::size_t>(nr) * static_cast<std::size_t>(W);

          auto scan_interval = [&](int a, int b) {
            std::size_t j = base + static_cast<std::size_t>(a);
            const std::size_t j_end = base + static_cast<std::size_t>(b);
            for (; j <= j_end; ++j) {
              if (!transit[j]) continue;
              if (visited_test_and_set(j, visited_bits, visited_word_stamp, cur_stamp)) continue;
              if (habitat[j]) g_sp += static_cast<double>(hab_weight[j]);
              next.push_back(j);
            }
          };

          for (std::size_t kk = 1; kk < ints.size(); ++kk) {
            const int a = ints[kk].first;
            const int b = ints[kk].second;
            if (a <= cur_b + 1) {
              if (b > cur_b) cur_b = b;
            } else {
              scan_interval(cur_a, cur_b);
              cur_a = a; cur_b = b;
            }
          }
          scan_interval(cur_a, cur_b);
        }
      }

      frontier.swap(next);
    }

    result[sp] = g_sp;
  } // sp
} // omp parallel

return result;
}
