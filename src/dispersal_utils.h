//------------------------------ Dispersal-aware G function --------------------

#ifndef DISPERSAL_UTILS_H
#define DISPERSAL_UTILS_H

#include <vector>
#include <cstdint>
#include "offsets.h"
#include "species_plan.h"

/* Function to calculate G accounting for dispersal through the landscape of habitat patches
 * - uses BFS approach for speed
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
    const SpeciesPlan& species_plan
);

#endif  // DISPERSAL_UTILS_H
