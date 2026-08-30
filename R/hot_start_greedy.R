#' Generate a Greedy Hot-Start Solution for EBREL
#'
#' Constructs a constrained initial action configuration for the EBREL
#' optimiser using a species-prioritised greedy allocation heuristic.
#'
#' Species are ranked according to their current range size and the number of
#' available habitat creation opportunities within a dispersal-based buffer.
#' A subset of the most constrained species is then selected using
#' `n_seed_species`, `max_opportunity_ratio`, or both. Habitat creation actions
#' are allocated sequentially for these seeded species.
#'
#' For each seeded species, the function identifies habitat creation
#' opportunities within a buffer defined as the species dispersal distance
#' multiplied by `buffer_multiplier`. Existing allocations that benefit the
#' species are retained, and additional cells are selected by proximity to the
#' current distribution. Where multiple suitable habitat actions are available,
#' the function favours the habitat type with the fewest available creation
#' opportunities across the landscape.
#'
#' Optionally, redundant creation allocations are pruned while preserving the
#' target contribution already achieved for the seeded species. The function
#' can then allocate habitat improvement actions within the existing
#' distributions of seeded species whose targets remain unmet.
#'
#' The returned raster uses the internal EBREL action coding: `-1` denotes no
#' action and `0:(n_actions - 1)` denote zero-indexed actions. Habitat creation
#' actions occupy `0:(n_habitats - 1)`, while the final action code,
#' `n_actions - 1`, represents habitat improvement.
#'
#' @param ebrel_sim_data A named list containing the EBREL input data. It must
#'   include:
#'   \describe{
#'     \item{`SD`}{A multilayer [terra::SpatRaster] containing one binary
#'       current-distribution layer per species. Occupied cells must be coded
#'       as `1`.}
#'     \item{`C`}{A multilayer [terra::SpatRaster] containing one cost layer
#'       per action. The final layer is assumed to represent habitat
#'       improvement; preceding layers represent habitat creation actions.}
#'     \item{`D`}{A numeric vector of species dispersal distances, with one
#'       value per species and in the same distance units as the raster CRS.}
#'     \item{`SxH`}{A species-by-habitat matrix identifying suitable habitat
#'       creation actions. Non-zero entries indicate habitat association.}
#'     \item{`O`}{A numeric vector containing one target value per species.
#'       Its interpretation is controlled by `target_mode`.}
#'   }
#'
#' @param sentinel Numeric scalar used in the cost raster to identify
#'   unavailable actions. Defaults to `1e10`.
#'
#' @param buffer_multiplier Numeric scalar multiplying each species dispersal
#'   distance when defining its focal allocation area. For example, a value of
#'   `2` searches within twice the supplied dispersal distance. Defaults to `2`.
#'
#' @param target_mode Character string specifying how values in `O` are
#'   interpreted. One of:
#'   \describe{
#'     \item{`"range_multiplier"`}{The requested target is calculated as
#'       `ceiling(current_range * O)`.}
#'     \item{`"absolute_cells"`}{Values in `O` are interpreted directly as the
#'       requested target contribution.}
#'   }
#'
#' @param rank_weights Named numeric vector giving the relative weights applied
#'   to species rarity and habitat opportunity when calculating the constraint
#'   score. It must contain elements named `"rarity"` and `"opportunity"`.
#'   Values are normalised internally to sum to one. Defaults to
#'   `c(rarity = 0.5, opportunity = 0.5)`.
#'
#' @param n_seed_species Optional non-negative integer giving the maximum number
#'   of species used to construct the hot start. Species are selected after
#'   ranking and after applying `max_opportunity_ratio`, where supplied.
#'   `NULL` applies no numerical cap. Defaults to `NULL`.
#'
#' @param max_opportunity_ratio Optional non-negative numeric threshold used to
#'   restrict the hot start to species with limited creation opportunities. The
#'   opportunity ratio is calculated as the number of available creation cells
#'   within the dispersal buffer divided by the requested target. Species with
#'   ratios less than or equal to this value are eligible for seeding.
#'   `NULL` applies no ratio filter. Defaults to `NULL`.
#'
#' @param allocate_improvement Logical scalar indicating whether habitat
#'   improvement actions should be allocated after creation allocation and
#'   pruning. Improvement is considered only for seeded species whose targets
#'   remain unmet. Defaults to `FALSE`.
#'
#' @param improvement_weight Positive numeric scalar giving the target
#'   contribution of one habitat improvement action. For example, a value of
#'   `0.5` means that two improvement cells contribute one target unit. This
#'   should normally match the improvement entry in `action_weight` supplied to
#'   [run_ebrel_R()]. Defaults to `0.1`.
#'
#' @param prune Logical scalar indicating whether redundant selected creation
#'   cells should be removed after greedy allocation. A cell is removed only
#'   where doing so does not reduce any benefiting seeded species below its
#'   preservation floor. Defaults to `TRUE`.
#'
#' @param seed Integer random seed used to resolve exact ranking, habitat
#'   selection, and pruning ties. Defaults to `123`.
#'
#' @return A single-layer [terra::SpatRaster] named `"X0"` containing the
#'   initial action configuration. Values use the internal EBREL coding:
#'   \describe{
#'     \item{`-1`}{No action selected.}
#'     \item{`0:(n_habitats - 1)`}{Zero-indexed habitat creation action.}
#'     \item{`n_actions - 1`}{Zero-indexed habitat improvement action.}
#'   }
#'
#' @details
#' Species are ranked using a weighted geometric mean of two rank-based scores:
#' current range rarity and scarcity of habitat opportunities within the
#' dispersal buffer. Species with small ranges and few available habitat
#' creation opportunities receive higher constraint scores.
#'
#' The opportunity ratio is calculated using the requested target rather than
#' the feasible target:
#'
#' \deqn{
#'   \mathrm{opportunity\ ratio} =
#'   \frac{\mathrm{available\ creation\ cells}}
#'        {\mathrm{requested\ target}}
#' }
#'
#' Values below `1` indicate that the requested target cannot be met using the
#' available creation opportunities. A value of `1` indicates exactly enough
#' opportunity, while larger values indicate increasing flexibility.
#'
#' When both `max_opportunity_ratio` and `n_seed_species` are supplied, the
#' ratio threshold is applied first. The function then retains at most
#' `n_seed_species` from the eligible species, preserving the constraint
#' ranking.
#'
#' For each seeded species, previously selected creation cells are counted where
#' their assigned habitat is suitable for that species and lies within its focal
#' buffer. Additional cells are selected from unallocated suitable cells,
#' ordered first by distance from the current distribution and then by
#' landscape-wide habitat rarity.
#'
#' Requested creation targets are capped at the number of available creation
#' opportunities when constructing the initial creation solution. Consequently,
#' the creation stage may retain a shortfall where the requested target is
#' infeasible.
#'
#' When `allocate_improvement = TRUE`, remaining shortfalls are evaluated after
#' creation allocation and pruning. Improvement cells are drawn from available,
#' currently unallocated cells within each seeded species' existing
#' distribution. Existing improvement allocations may contribute to multiple
#' seeded species where their distributions overlap.
#'
#' Improvement cells are allocated according to their costs in the final layer
#' of `C`. The number of cells required is determined from the remaining target
#' contribution and `improvement_weight`.
#'
#' The function is intended to provide a sparse, ecologically plausible
#' starting configuration rather than a complete solution for all species. The
#' simulated-annealing optimiser remains responsible for resolving remaining
#' targets, costs, and spatial configuration trade-offs.
#'
#' @seealso [run_ebrel_R()]
#'
#' @examples
#' \dontrun{
#' # Seed at most ten species with no more than twice the creation
#' # opportunity required to meet their requested target.
#' X0 <- hot_start_greedy(
#'   ebrel_sim_data,
#'   sentinel = 1e10,
#'   buffer_multiplier = 2,
#'   target_mode = "range_multiplier",
#'   rank_weights = c(
#'     rarity = 0.5,
#'     opportunity = 0.5
#'   ),
#'   n_seed_species = 10,
#'   max_opportunity_ratio = 2,
#'   allocate_improvement = TRUE,
#'   improvement_weight = 0.1,
#'   prune = TRUE,
#'   seed = 123
#' )
#'
#' terra::plot(X0)
#'
#' result <- run_ebrel_R(
#'   E = ebrel_sim_data$E,
#'   C = ebrel_sim_data$C,
#'   SD = ebrel_sim_data$SD,
#'   D = ebrel_sim_data$D,
#'   SxH = ebrel_sim_data$SxH,
#'   O = ebrel_sim_data$O,
#'   action_weight = c(
#'     rep(1, terra::nlyr(ebrel_sim_data$C) - 1L),
#'     0.5
#'   ),
#'   LM = ebrel_sim_data$LM,
#'   X0 = X0
#' )
#' }
#'
#' @export
hot_start_greedy <- function(
    ebrel_sim_data,
    sentinel = 1e10,
    buffer_multiplier = 2,
    target_mode = c("range_multiplier", "absolute_cells"),
    rank_weights = c(rarity = 0.5, opportunity = 0.5),
    n_seed_species = NULL,
    max_opportunity_ratio = NULL,
    allocate_improvement = FALSE,
    improvement_weight = 0.1,
    prune = TRUE,
    seed = 123
) {
  
  stopifnot(
    inherits(ebrel_sim_data$SD, "SpatRaster"),
    inherits(ebrel_sim_data$C, "SpatRaster")
  )
  
  target_mode <- match.arg(target_mode)
  set.seed(seed)
  
  # ---- Inputs ----
  
  Dmax <- ebrel_sim_data$D
  SD   <- ebrel_sim_data$SD
  C    <- ebrel_sim_data$C
  SxH  <- as.matrix(ebrel_sim_data$SxH)
  O    <- ebrel_sim_data$O
  
  n_species  <- terra::nlyr(SD)
  n_actions  <- terra::nlyr(C)
  n_habitats <- n_actions - 1L
  n_cells    <- terra::ncell(C)
  
  stopifnot(
    terra::compareGeom(SD, C, stopOnError = FALSE),
    length(Dmax) == n_species,
    length(O) == n_species,
    nrow(SxH) == n_species,
    ncol(SxH) == n_habitats
  )
  
  if (is.null(names(C))) {
    names(C) <- paste0("C", seq_len(n_habitats))
  }
  
  if (is.null(names(SD))) {
    names(SD) <- paste0("species_", seq_len(n_species))
  }
  
  rank_weights <- rank_weights / sum(rank_weights)
  
  
  # ---- Habitat availability and habitat rarity ----
  
  C_values <- terra::values(C, mat = TRUE)
  
  # The last C layer is habitat improvement. The hot-start allocation below
  # currently assigns creation actions only, so retain only creation layers.
  creation_costs <- C_values[
    ,
    seq_len(n_habitats),
    drop = FALSE
  ]
  
  # Match this tolerance to the C++ approx_equal() implementation.
  tolerance <- sqrt(.Machine$double.eps)
  
  comparison_scale <- pmax(
    abs(creation_costs),
    abs(sentinel),
    1
  )
  
  creation_unavailable <- (
    abs(creation_costs - sentinel) <=
      tolerance * comparison_scale
  )
  
  habitat_available <- !creation_unavailable
  
  # The final action layer represents improvement of existing habitat.
  improvement_costs <- C_values[, n_actions]
  improvement_scale <- pmax(abs(improvement_costs), abs(sentinel), 1)
  improvement_unavailable <- (
    abs(improvement_costs - sentinel) <=
      tolerance * improvement_scale
  )
  improvement_available <- !improvement_unavailable
  
  if (
    !is.logical(allocate_improvement) ||
    length(allocate_improvement) != 1L ||
    is.na(allocate_improvement)
  ) {
    stop("`allocate_improvement` must be TRUE or FALSE.")
  }
  
  if (
    !is.numeric(improvement_weight) ||
    length(improvement_weight) != 1L ||
    is.na(improvement_weight) ||
    !is.finite(improvement_weight) ||
    improvement_weight <= 0
  ) {
    stop("`improvement_weight` must be a single positive numeric value.")
  }
  
  # Habitat rarity is measured across available creation opportunities.
  habitat_extent <- colSums(habitat_available)
  
  habitat_info <- data.frame(
    habitat_id = seq_len(n_habitats),
    habitat = names(C)[seq_len(n_habitats)],
    available_cells = habitat_extent,
    stringsAsFactors = FALSE
  )
  
  # ---- Derive species-specific spatial information ----
  
  species_info <- data.frame(
    species_id = seq_len(n_species),
    species_name = names(SD),
    current_range = integer(n_species),
    habitat_opp = integer(n_species),
    improvement_opp = integer(n_species),
    target_cells = integer(n_species),
    n_associated_habitats = integer(n_species),
    stringsAsFactors = FALSE
  )
  
  # Retain information needed during allocation.
  species_distance <- vector("list", n_species)
  species_focal    <- vector("list", n_species)
  species_habitats <- vector("list", n_species)
  species_occupied <- vector("list", n_species)
  
  for (i in seq_len(n_species)) {
    
    sd_values <- terra::values(SD[[i]], mat = FALSE)
    
    occupied <- which(
      !is.na(sd_values) &
        sd_values == 1
    )
    
    species_info$current_range[i] <- length(occupied)
    species_occupied[[i]] <- occupied
    
    associated_habitats <- which(SxH[i, ] == 1)
    
    species_habitats[[i]] <- associated_habitats
    
    species_info$n_associated_habitats[i] <-
      length(associated_habitats)
    
    # Improvement opportunities occur only within the current distribution.
    species_info$improvement_opp[i] <- sum(
      improvement_available[occupied]
    )
    
    # Species without occupied cells cannot seed either creation or
    # improvement actions.
    if (length(occupied) == 0L) {
      species_distance[[i]] <- rep(NA_real_, n_cells)
      species_focal[[i]] <- rep(FALSE, n_cells)
      species_info$habitat_opp[i] <- 0L
      next
    }
    
    # Raster with only occupied cells retained.
    occupied_rast <- terra::rast(SD[[i]])
    terra::values(occupied_rast) <- NA_real_
    occupied_rast[occupied] <- 1
    
    # Euclidean distance to the nearest currently occupied cell.
    distance_rast <- terra::distance(occupied_rast)
    distance_values <- terra::values(distance_rast, mat = FALSE)
    
    species_distance[[i]] <- distance_values
    
    buffer_width <- Dmax[i] * buffer_multiplier
    
    focal_cells <- (
      !is.na(distance_values) &
        distance_values <= buffer_width
    )
    
    species_focal[[i]] <- focal_cells
    
    # Cells containing at least one suitable habitat-creation opportunity.
    if (length(associated_habitats) > 0L) {
      opportunity_cells <- focal_cells &
        rowSums(
          habitat_available[, associated_habitats, drop = FALSE]
        ) > 0
      species_info$habitat_opp[i] <- sum(opportunity_cells)
    } else {
      species_info$habitat_opp[i] <- 0L
    }
  }
  
  # ---- Define species targets ----
  
  if (target_mode == "range_multiplier") {
    
    # O is interpreted as a target number of action cells per currently
    # occupied cell, e.g. O = 0.5 gives ceiling(current_range * 0.5).
    species_info$target_cells <- ceiling(
      species_info$current_range * O
    )
    
  } else {
    
    # O is interpreted directly as a number of cells.
    species_info$target_cells <- ceiling(O)
  }
  
  species_info$target_cells[
    species_info$current_range == 0
  ] <- 0L
  
  # Targets cannot exceed the available opportunity.
  species_info$feasible_target <- pmin(
    species_info$target_cells,
    species_info$habitat_opp
  )
  
  species_info$target_feasible <- (
    species_info$habitat_opp >= species_info$target_cells
  )
  
  # Opportunity relative to the requested target. Low values identify species
  # with relatively little flexibility in where their target can be met.
  # Species with no positive target receive Inf and are not selected by a
  # finite max_opportunity_ratio threshold.
  species_info$opportunity_ratio <- Inf
  positive_target <- species_info$target_cells > 0L
  species_info$opportunity_ratio[positive_target] <-
    species_info$habitat_opp[positive_target] /
    species_info$target_cells[positive_target]
  
  # ---- Balanced constraint ranking ----
  
  # A rank-based score avoids sensitivity to very large differences in the
  # raw range-size and opportunity variables.
  rarity_rank <- rank(
    species_info$current_range,
    ties.method = "average"
  )
  
  opportunity_rank <- rank(
    species_info$habitat_opp,
    ties.method = "average"
  )
  
  # Small range and small opportunity receive high scores.
  species_info$rarity_score <-
    (n_species - rarity_rank + 1) / n_species
  
  species_info$opportunity_score <-
    (n_species - opportunity_rank + 1) / n_species
  
  # Weighted geometric mean requires both constraint dimensions to matter.
  species_info$constraint_score <-
    species_info$rarity_score^rank_weights["rarity"] *
    species_info$opportunity_score^rank_weights["opportunity"]
  
  # Random term only resolves exact ranking ties.
  species_info$tie_break <- runif(n_species)
  
  species_info <- species_info[
    order(
      -species_info$constraint_score,
      species_info$current_range,
      species_info$habitat_opp,
      species_info$tie_break
    ),
  ]
  
  species_info$allocation_order <- seq_len(n_species)
  
  # ---- Select species used to construct the hot start ----
  
  # Only species with a positive target and at least one available creation
  # or improvement opportunity can contribute to the initial solution.
  seed_eligible <- (
    species_info$target_cells > 0L &
      (
        species_info$feasible_target > 0L |
          (
            allocate_improvement &
              species_info$improvement_opp > 0L
          )
      )
  )
  
  # Optionally retain only species with relatively few habitat opportunities
  # per requested target cell.
  if (!is.null(max_opportunity_ratio)) {
    
    if (
      !is.numeric(max_opportunity_ratio) ||
      length(max_opportunity_ratio) != 1L ||
      is.na(max_opportunity_ratio) ||
      !is.finite(max_opportunity_ratio) ||
      max_opportunity_ratio < 0
    ) {
      stop(
        "`max_opportunity_ratio` must be NULL or a single ",
        "non-negative finite numeric value."
      )
    }
    
    seed_eligible <- (
      seed_eligible &
        species_info$opportunity_ratio <= max_opportunity_ratio
    )
  }
  
  # species_info is already ordered from most to least constrained.
  seed_species_info <- species_info[
    seed_eligible,
    ,
    drop = FALSE
  ]
  
  # Optionally cap the number of seeded species after applying the opportunity
  # ratio filter.
  if (!is.null(n_seed_species)) {
    
    if (
      !is.numeric(n_seed_species) ||
      length(n_seed_species) != 1L ||
      is.na(n_seed_species) ||
      !is.finite(n_seed_species) ||
      n_seed_species < 0 ||
      n_seed_species != as.integer(n_seed_species)
    ) {
      stop(
        "`n_seed_species` must be NULL or a single non-negative integer."
      )
    }
    
    n_seed_species <- min(
      as.integer(n_seed_species),
      nrow(seed_species_info)
    )
    
    seed_species_info <- seed_species_info[
      seq_len(n_seed_species),
      ,
      drop = FALSE
    ]
  }
  
  seed_species_ids <- seed_species_info$species_id
  
  # ---- Greedy allocation ----
  
  # Each element is either NA or a global habitat-layer index.
  selected_habitat <- rep(NA_integer_, n_cells)
  
  # Improvement is stored separately because selected_habitat contains
  # creation-habitat ids rather than general action-layer ids.
  selected_improvement <- rep(FALSE, n_cells)
  
  allocation_records <- vector("list", nrow(seed_species_info))
  
  for (j in seq_len(nrow(seed_species_info))) {
    
    i <- seed_species_info$species_id[j]
    
    associated_habitats <- species_habitats[[i]]
    focal_cells <- species_focal[[i]]
    distance_values <- species_distance[[i]]
    
    target <- seed_species_info$feasible_target[j]
    
    if (
      target <= 0 ||
      length(associated_habitats) == 0 ||
      !any(focal_cells)
    ) {
      
      allocation_records[[j]] <- data.frame(
        species_id = i,
        species_name = names(SD)[i],
        allocation_order = j,
        target_cells = target,
        previously_contributing = 0L,
        newly_allocated = 0L,
        achieved_cells = 0L,
        shortfall = target,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    # Existing allocations can contribute where:
    # 1. the cell lies within the species focal area; and
    # 2. its allocated habitat is associated with the species.
    existing_contributing <- which(
      focal_cells &
        !is.na(selected_habitat) &
        selected_habitat %in% associated_habitats
    )
    
    n_existing <- length(existing_contributing)
    n_required <- max(0L, target - n_existing)
    
    selected_now <- integer(0)
    
    if (n_required > 0) {
      
      # Candidate cells must be unallocated and contain at least one suitable
      # habitat opportunity for this species.
      candidate_cells <- which(
        focal_cells &
          is.na(selected_habitat) &
          rowSums(
            habitat_available[
              ,
              associated_habitats,
              drop = FALSE
            ]
          ) > 0
      )
      
      if (length(candidate_cells) > 0) {
        
        candidate_habitat <- integer(length(candidate_cells))
        candidate_habitat_extent <- numeric(length(candidate_cells))
        
        for (k in seq_along(candidate_cells)) {
          
          cell_id <- candidate_cells[k]
          
          available_here <- associated_habitats[
            habitat_available[
              cell_id,
              associated_habitats
            ]
          ]
          
          rarity_here <- habitat_extent[available_here]
          minimum_extent <- min(rarity_here)
          
          # Habitat types with the smallest landscape-wide extent.
          rarest_here <- available_here[
            rarity_here == minimum_extent
          ]
          
          # Random selection only where habitat rarity is tied.
          candidate_habitat[k] <- rarest_here[
            sample.int(length(rarest_here), 1L)
          ]
          candidate_habitat_extent[k] <- minimum_extent
        }
        
        candidate_table <- data.frame(
          cell = candidate_cells,
          habitat_id = candidate_habitat,
          distance_m = distance_values[candidate_cells],
          habitat_extent = candidate_habitat_extent,
          tie_break = runif(length(candidate_cells))
        )
        
        # Distance is the primary criterion. Habitat rarity and a random
        # number resolve cells at equal or very similar distance.
        candidate_table <- candidate_table[
          order(
            candidate_table$distance_m,
            candidate_table$habitat_extent,
            candidate_table$tie_break
          ),
        ]
        
        n_take <- min(n_required, nrow(candidate_table))
        
        selected_now <- candidate_table$cell[
          seq_len(n_take)
        ]
        
        selected_habitat[selected_now] <-
          candidate_table$habitat_id[
            seq_len(n_take)
          ]
      }
    }
    
    achieved <- n_existing + length(selected_now)
    
    allocation_records[[j]] <- data.frame(
      species_id = i,
      species_name = names(SD)[i],
      allocation_order = j,
      target_cells = target,
      previously_contributing = n_existing,
      newly_allocated = length(selected_now),
      achieved_cells = achieved,
      shortfall = max(0L, target - achieved),
      stringsAsFactors = FALSE
    )
  }
  
  allocation_summary <- if (length(allocation_records) > 0L) {
    do.call(rbind, allocation_records)
  } else {
    data.frame(
      species_id = integer(0),
      species_name = character(0),
      allocation_order = integer(0),
      target_cells = integer(0),
      previously_contributing = integer(0),
      newly_allocated = integer(0),
      achieved_cells = integer(0),
      shortfall = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Prune allocations that are redundant across all species ----
  
  n_pruned <- 0L
  
  if (
    prune &&
    length(seed_species_ids) > 0L &&
    any(!is.na(selected_habitat))
  ) {
    
    # Identify species benefiting from each selected action.
    contribution_matrix <- matrix(
      FALSE,
      nrow = n_cells,
      ncol = length(seed_species_ids)
    )
    
    for (k in seq_along(seed_species_ids)) {
      
      i <- seed_species_ids[k]
      associated_habitats <- species_habitats[[i]]
      
      if (length(associated_habitats) == 0L) {
        next
      }
      
      contribution_matrix[, k] <- (
        species_focal[[i]] &
          !is.na(selected_habitat) &
          selected_habitat %in% associated_habitats
      )
    }
    
    achieved_before_pruning <- colSums(contribution_matrix)
    
    requested_targets <- seed_species_info$feasible_target
    
    # For infeasible or unmet species, preserve the achieved amount rather
    # than allowing pruning to worsen the shortfall.
    preservation_floor <- pmin(
      requested_targets,
      achieved_before_pruning
    )
    
    current_achievement <- achieved_before_pruning
    
    selected_cells <- which(!is.na(selected_habitat))
    
    # Test cells contributing to the fewest species first.
    cell_contribution_count <- rowSums(
      contribution_matrix[selected_cells, , drop = FALSE]
    )
    
    prune_order <- selected_cells[
      order(
        cell_contribution_count,
        runif(length(selected_cells))
      )
    ]
    
    for (cell_id in prune_order) {
      
      benefiting_species <- which(
        contribution_matrix[cell_id, ]
      )
      
      if (length(benefiting_species) == 0) {
        
        selected_habitat[cell_id] <- NA_integer_
        n_pruned <- n_pruned + 1L
        next
      }
      
      can_remove <- all(
        current_achievement[benefiting_species] - 1 >=
          preservation_floor[benefiting_species]
      )
      
      if (can_remove) {
        
        selected_habitat[cell_id] <- NA_integer_
        
        current_achievement[benefiting_species] <-
          current_achievement[benefiting_species] - 1L
        
        contribution_matrix[cell_id, ] <- FALSE
        n_pruned <- n_pruned + 1L
      }
    }
  }
  
  # ---- Optionally allocate improvement to occupied cells ----
  
  improvement_records <- vector("list", nrow(seed_species_info))
  
  if (allocate_improvement && length(seed_species_ids) > 0L) {
    
    # Improvement is allocated only after creation has been pruned. Seeded
    # species are processed in the same constrained-species order used above.
    # An improvement action may contribute to more than one species where their
    # existing distributions overlap.
    for (j in seq_len(nrow(seed_species_info))) {
      
      i <- seed_species_info$species_id[j]
      associated_habitats <- species_habitats[[i]]
      occupied <- species_occupied[[i]]
      
      creation_contribution <- sum(
        species_focal[[i]] &
          !is.na(selected_habitat) &
          selected_habitat %in% associated_habitats
      )
      
      existing_improvement_cells <- occupied[
        selected_improvement[occupied]
      ]
      
      improvement_contribution <-
        length(existing_improvement_cells) * improvement_weight
      
      requested_target <- seed_species_info$target_cells[j]
      remaining_target <- max(
        0,
        requested_target -
          creation_contribution -
          improvement_contribution
      )
      
      selected_now <- integer(0)
      
      if (remaining_target > 0 && length(occupied) > 0L) {
        
        # A cell can contain only one action. Candidate improvement cells must
        # therefore be occupied by the species, currently unallocated, not
        # already selected for improvement, and available in the final C layer.
        candidate_cells <- occupied[
          is.na(selected_habitat[occupied]) &
            !selected_improvement[occupied] &
            improvement_available[occupied]
        ]
        
        if (length(candidate_cells) > 0L) {
          
          # Prefer lower-cost improvement cells; random values resolve exact
          # cost ties while retaining reproducibility through seed.
          candidate_table <- data.frame(
            cell = candidate_cells,
            cost = improvement_costs[candidate_cells],
            tie_break = runif(length(candidate_cells))
          )
          
          candidate_table <- candidate_table[
            order(
              candidate_table$cost,
              candidate_table$tie_break
            ),
            ,
            drop = FALSE
          ]
          
          n_required <- ceiling(
            remaining_target / improvement_weight
          )
          n_take <- min(n_required, nrow(candidate_table))
          
          selected_now <- candidate_table$cell[
            seq_len(n_take)
          ]
          selected_improvement[selected_now] <- TRUE
        }
      }
      
      final_improvement_contribution <- (
        sum(selected_improvement[occupied]) * improvement_weight
      )
      final_creation_contribution <- sum(
        species_focal[[i]] &
          !is.na(selected_habitat) &
          selected_habitat %in% associated_habitats
      )
      final_contribution <-
        final_creation_contribution + final_improvement_contribution
      
      improvement_records[[j]] <- data.frame(
        species_id = i,
        species_name = names(SD)[i],
        allocation_order = j,
        target_cells = requested_target,
        creation_contribution = final_creation_contribution,
        previous_improvement_contribution = improvement_contribution,
        newly_improved_cells = length(selected_now),
        improvement_contribution = final_improvement_contribution,
        achieved_contribution = final_contribution,
        shortfall = max(0, requested_target - final_contribution),
        stringsAsFactors = FALSE
      )
    }
  }
  
  improvement_summary <- if (
    allocate_improvement &&
    length(improvement_records) > 0L
  ) {
    do.call(rbind, improvement_records)
  } else {
    data.frame(
      species_id = integer(0),
      species_name = character(0),
      allocation_order = integer(0),
      target_cells = numeric(0),
      creation_contribution = numeric(0),
      previous_improvement_contribution = numeric(0),
      newly_improved_cells = integer(0),
      improvement_contribution = numeric(0),
      achieved_contribution = numeric(0),
      shortfall = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Final species outcomes ----
  
  final_achievement <- numeric(n_species)
  
  for (i in seq_len(n_species)) {
    
    associated_habitats <- species_habitats[[i]]
    
    creation_contribution <- if (length(associated_habitats) > 0L) {
      sum(
        species_focal[[i]] &
          !is.na(selected_habitat) &
          selected_habitat %in% associated_habitats
      )
    } else {
      0
    }
    
    improvement_contribution <- (
      sum(selected_improvement[species_occupied[[i]]]) *
        improvement_weight
    )
    
    final_achievement[i] <-
      creation_contribution + improvement_contribution
  }
  
  final_target <- species_info$target_cells[
    match(
      seq_len(n_species),
      species_info$species_id
    )
  ]
  
  species_outcomes <- data.frame(
    species_id = seq_len(n_species),
    species_name = names(SD),
    target_cells = final_target,
    achieved_cells = final_achievement,
    shortfall = pmax(0L, final_target - final_achievement),
    target_met = final_achievement >= final_target,
    stringsAsFactors = FALSE
  )
  
  species_outcomes <- merge(
    species_info,
    species_outcomes,
    by = c(
      "species_id",
      "species_name",
      "target_cells"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  species_outcomes <- species_outcomes[
    order(species_outcomes$allocation_order),
  ]
  
  # ---- Validate final action assignments against C ----
  
  X0_values <- rep.int(-1L, n_cells)
  
  creation_cells <- which(!is.na(selected_habitat))
  improvement_cells <- which(selected_improvement)
  
  if (length(creation_cells) > 0L) {
    X0_values[creation_cells] <-
      selected_habitat[creation_cells] - 1L
    
    selected_layers <- X0_values[creation_cells] + 1L
    
    selected_unavailable <- creation_unavailable[
      cbind(creation_cells, selected_layers)
    ]
    
    selected_costs <- creation_costs[
      cbind(creation_cells, selected_layers)
    ]
    
    if (any(selected_unavailable)) {
      bad <- which(selected_unavailable)[1L]
      
      stop(
        "Hot start selected an unavailable creation action at cell ",
        creation_cells[bad],
        ", action ",
        X0_values[creation_cells[bad]],
        ", habitat layer ",
        selected_layers[bad],
        ", cost ",
        format(selected_costs[bad], digits = 17)
      )
    }
  }
  
  if (length(improvement_cells) > 0L) {
    
    if (any(!is.na(selected_habitat[improvement_cells]))) {
      stop(
        "Hot start assigned both creation and improvement to the same cell."
      )
    }
    
    if (any(improvement_unavailable[improvement_cells])) {
      bad <- which(improvement_unavailable[improvement_cells])[1L]
      stop(
        "Hot start selected an unavailable improvement action at cell ",
        improvement_cells[bad],
        ", action ",
        n_actions - 1L,
        ", cost ",
        format(improvement_costs[improvement_cells[bad]], digits = 17)
      )
    }
    
    # The final C layer is the improvement action. Convert its one-based R
    # layer index to the zero-based internal EBREL action code.
    X0_values[improvement_cells] <- n_actions - 1L
  }
  
  if (any(X0_values < -1L | X0_values >= n_actions)) {
    stop(
      "Hot start produced action codes outside -1:",
      n_actions - 1L
    )
  }
  
  # ---- Output rasters ----
  
  action_raster <- terra::rast(C[[1]])
  terra::values(action_raster) <- X0_values
  names(action_raster) <- "X0"
  
  # ---- Return ----
  
  action_raster
  
}
