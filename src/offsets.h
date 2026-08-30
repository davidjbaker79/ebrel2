#pragma once

#include <vector>
#include <utility>
#include <cstdint>

// Use int16_t to reduce memory/bandwidth
using Offset  = std::pair<int16_t,int16_t>;
using Offsets = std::vector<Offset>;
using OffsetsCache = std::vector<Offsets>;

// Build offsets for one dispersal distance d (Euclidean disk)
Offsets make_offsets_round(int d);

// Build cache for all d = 0..universal_disp_thres
OffsetsCache build_offsets_cache_round(int universal_disp_thres);

// For each dr, store the inclusive dc range [dc_min, dc_max]
struct RowRun { int16_t dr; int16_t dc_min; int16_t dc_max; };
using RowRuns = std::vector<RowRun>;
using RowRunsCache = std::vector<RowRuns>;

//  row runs grouping from offsets
RowRuns group_offsets_by_dr_to_runs(const Offsets& off);

// Build row runs cache
RowRunsCache build_rowruns_cache_from_offsets(const OffsetsCache& offsets_cache);