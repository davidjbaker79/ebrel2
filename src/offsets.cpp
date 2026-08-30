#include "offsets.h"

#include <cmath>      // std::sqrt, std::floor
#include <cstddef>    // std::size_t
#include <algorithm>  // std::max

Offsets make_offsets_round(int d)
{
  Offsets off;
  if (d <= 0) return off;
  
  // Rough reserve: area of disk ~ pi*d^2, plus a little slack
  off.reserve(static_cast<std::size_t>(3.2 * d * d));
  
  const int d2 = d * d;
  
  for (int dr = -d; dr <= d; ++dr) {
    const int dr2 = dr * dr;
    for (int dc = -d; dc <= d; ++dc) {
      if (dr == 0 && dc == 0) continue;
      if (dr2 + dc * dc <= d2) {
        off.emplace_back(static_cast<int16_t>(dr),
                         static_cast<int16_t>(dc));
      }
    }
  }
  
  return off;
}

// Build offsets cache
OffsetsCache build_offsets_cache_round(int universal_disp_thres)
{
  if (universal_disp_thres < 0) universal_disp_thres = 0;
  
  OffsetsCache cache(static_cast<std::size_t>(universal_disp_thres) + 1);
  
  for (int d = 0; d <= universal_disp_thres; ++d) {
    cache[static_cast<std::size_t>(d)] = make_offsets_round(d);
  }
  
  return cache;
}

// Row runs grouping from offsets
RowRuns group_offsets_by_dr_to_runs(const Offsets& off) {
  RowRuns runs;
  if (off.empty()) return runs;
  // dr -> (min_dc, max_dc)
  // For a disk, for each dr there will be a contiguous range.
  // Use a small map keyed by dr. Since dr is in [-d, d], use vector indexing.
  int16_t min_dr = 32767, max_dr = -32768;
  for (auto &p : off) { min_dr = std::min(min_dr, p.first); max_dr = std::max(max_dr, p.first); }
  
  const int span = (max_dr - min_dr + 1);
  std::vector<int16_t> lo(span,  32767);
  std::vector<int16_t> hi(span, -32768);
  std::vector<uint8_t> has(span, 0u);
  for (auto &p : off) {
    const int idx = p.first - min_dr;
    has[idx] = 1u;
    lo[idx] = std::min(lo[idx], p.second);
    hi[idx] = std::max(hi[idx], p.second);
  }
  runs.reserve(span);
  for (int i = 0; i < span; ++i) {
    if (!has[i]) continue;
    runs.push_back(RowRun{ static_cast<int16_t>(min_dr + i), lo[i], hi[i] });
  }
  return runs;
}

// Build row runs cache
RowRunsCache build_rowruns_cache_from_offsets(const OffsetsCache& offsets_cache) {
  RowRunsCache out(offsets_cache.size());
  for (std::size_t d = 0; d < offsets_cache.size(); ++d) {
    out[d] = group_offsets_by_dr_to_runs(offsets_cache[d]);
  }
  return out;
}


