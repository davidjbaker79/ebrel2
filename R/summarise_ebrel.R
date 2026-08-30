#' Summarise an EBREL optimisation result
#'
#' Creates a structured, user-friendly summary of an EBREL optimisation result.
#' The summary includes species-level target attainment, the relative
#' contribution of habitat creation and improvement actions, selected action
#' counts, objective weights, calibration diagnostics, and optimisation-run
#' statistics.
#'
#' Species-level target attainment is derived from \code{g_best},
#' \code{g_create_best}, and \code{g_improve_best}. By convention,
#' \code{g_best} represents remaining target shortfall, whereas
#' \code{g_create_best} and \code{g_improve_best} represent progress towards
#' the target attributable to habitat creation and habitat improvement,
#' respectively.
#'
#' @param x A list returned by \code{run_ebrel_cpp()} containing, at minimum,
#'   \code{H_best}, \code{iterations_run}, \code{accepted}, \code{proposals},
#'   \code{overall_acc}, \code{g_best}, \code{g_create_best}, and
#'   \code{g_improve_best}. If \code{x$X_best} is a \code{terra::SpatRaster},
#'   the summary also includes the number and proportion of selected cells for
#'   each action.
#'
#' @param species_names Optional character vector giving species names in the
#'   same order as \code{g_best}, \code{g_create_best}, and
#'   \code{g_improve_best}. If \code{NULL}, generic names such as
#'   \code{"Species 1"} are used.
#'
#' @param action_names Optional character vector giving labels for action IDs.
#'   The vector must be indexed by action ID plus one, because EBREL action IDs
#'   are zero-based. For example, action IDs \code{0}, \code{1}, and \code{2}
#'   would correspond to:
#'   \preformatted{
#'   c("Create habitat 1", "Create habitat 2", "Improve habitat")
#'   }
#'   Missing labels are replaced with generic labels such as \code{"Action 2"}.
#'
#' @param shortfall_tolerance Numeric value defining when a species target is
#'   considered met. Species with remaining target shortfall less than or equal
#'   to this value are classified as having met their target. Defaults to
#'   \code{1e-8}.
#'
#' @param print Logical. Should a readable summary be printed before the object
#'   is returned? Defaults to \code{TRUE}.
#'
#' @details
#' The species summary reports:
#'
#' \itemize{
#'   \item \code{creation_target_contribution}: target progress attributed to
#'   creation actions;
#'   \item \code{improvement_target_contribution}: target progress attributed
#'   to improvement actions;
#'   \item \code{total_target_attainment}: combined progress towards the target,
#'   capped at one;
#'   \item \code{remaining_target_shortfall}: unmet proportion of the target;
#'   \item \code{target_met}: whether the remaining shortfall is within
#'   \code{shortfall_tolerance};
#'   \item \code{target_balance_check}: the difference between total attainment
#'   and \code{1 - g_best}. This is intended as a numerical consistency check.
#' }
#'
#' Where available, the function also extracts objective calibration information
#' from \code{x$objective_calibration}, including component scaling values and
#' action-specific calibration diagnostics. Objective preference weights are
#' extracted from \code{x$objective_weights}.
#'
#' @return
#' A list of class \code{"ebrel_summary"} with the following components:
#'
#' \describe{
#'   \item{headline}{A named list of headline statistics, including the number
#'   and proportion of species meeting their targets, mean and median target
#'   attainment, total contributions from creation and improvement, and the
#'   number of selected cells where available.}
#'
#'   \item{run}{A one-row data frame summarising the best objective value,
#'   iteration count, proposals, accepted moves, acceptance rate, early stopping
#'   iteration, and iteration timing information where available.}
#'
#'   \item{objective}{A one-row data frame containing the objective preference
#'   weights \code{alpha}, \code{beta}, and \code{gamma}, if supplied in the
#'   EBREL output.}
#'
#'   \item{calibration}{A one-row data frame containing objective calibration
#'   scaling values and calibration sampling diagnostics, if available.}
#'
#'   \item{calibration_by_action}{An action-level calibration table, if
#'   available. Action labels are added when \code{action_names} is supplied.}
#'
#'   \item{actions}{A data frame reporting selected cell counts and the
#'   percentage of all selected cells associated with each action. This is
#'   \code{NULL} when \code{x$X_best} is unavailable or is not a
#'   \code{terra::SpatRaster}.}
#'
#'   \item{species}{A data frame containing species-level target attainment,
#'   shortfall, action-specific contributions, and target status.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run an EBREL optimisation
#' fit <- run_ebrel_cpp(...)
#'
#' # Create a readable summary
#' summary_fit <- summarise_ebrel(
#'   fit,
#'   species_names = species_names,
#'   action_names = c(
#'     "Create broadleaved woodland",
#'     "Create species-rich grassland",
#'     "Improve existing habitat"
#'   )
#' )
#'
#' # Inspect species that did not meet their targets
#' subset(summary_fit$species, !target_met)
#'
#' # Inspect the selected action mix
#' summary_fit$actions
#' }
#'
#' @export
#'
#' @seealso
#' \code{\link{print.ebrel_summary}} for a concise printed representation of
#' the returned summary object.
summarise_ebrel <- function(
    x,
    species_names = NULL,
    action_names = NULL,
    shortfall_tolerance = 1e-8,
    print = TRUE
) {
  stopifnot(is.list(x))
  
  required <- c(
    "H_best",
    "iterations_run",
    "accepted",
    "proposals",
    "overall_acc",
    "g_best",
    "g_create_best",
    "g_improve_best"
  )
  
  missing_fields <- setdiff(required, names(x))
  
  if (length(missing_fields) > 0) {
    stop(
      "EBREL output is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Species-level target outcomesa
  n_species <- length(x$g_best)
  
  if (is.null(species_names)) {
    species_names <- paste0("Species ", seq_len(n_species))
  }
  
  if (length(species_names) != n_species) {
    stop(
      "species_names must have the same length as g_best.",
      call. = FALSE
    )
  }
  
  target_shortfall <- pmax(0, pmin(1, as.numeric(x$g_best)))
  
  creation_progress <- pmax(
    0,
    as.numeric(x$g_create_best)
  )
  
  improvement_progress <- pmax(
    0,
    as.numeric(x$g_improve_best)
  )
  
  total_target_attainment <- pmin(
    1,
    creation_progress + improvement_progress
  )
  
  # This is useful as a consistency check. It should normally be approximately
  # zero, apart from numerical rounding or any future changes to g logic.
  attainment_difference <-
    total_target_attainment - (1 - target_shortfall)
  
  species_summary <- data.frame(
    species = species_names,
    
    creation_target_contribution =
      creation_progress,
    
    improvement_target_contribution =
      improvement_progress,
    
    total_target_attainment =
      total_target_attainment,
    
    target_attainment_percent =
      100 * total_target_attainment,
    
    remaining_target_shortfall =
      target_shortfall,
    
    target_shortfall_percent =
      100 * target_shortfall,
    
    target_met =
      target_shortfall <= shortfall_tolerance,
    
    target_balance_check =
      attainment_difference,
    
    stringsAsFactors = FALSE
  )
  
  # Action selection summary
  
  action_summary <- NULL
  
  if (!is.null(x$X_best) &&
      inherits(x$X_best, "SpatRaster")) {
    
    selected_actions <- terra::values(x$X_best, mat = FALSE)
    
    selected_actions <- selected_actions[
      !is.na(selected_actions) & selected_actions >= 0
    ]
    
    if (length(selected_actions) == 0) {
      
      action_summary <- data.frame(
        action_id = integer(0),
        action = character(0),
        selected_cells = integer(0),
        selected_percent_of_all_selected = numeric(0),
        stringsAsFactors = FALSE
      )
      
    } else {
      
      action_ids <- sort(unique(as.integer(selected_actions)))
      
      selected_counts <- vapply(
        action_ids,
        function(action_id) {
          sum(selected_actions == action_id)
        },
        numeric(1)
      )
      
      if (is.null(action_names)) {
        action_labels <- paste0("Action ", action_ids)
      } else {
        action_labels <- ifelse(
          action_ids + 1L <= length(action_names),
          action_names[action_ids + 1L],
          paste0("Action ", action_ids)
        )
      }
      
      action_summary <- data.frame(
        action_id = action_ids,
        action = action_labels,
        selected_cells = as.integer(selected_counts),
        selected_percent_of_all_selected =
          100 * selected_counts / sum(selected_counts),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Calibration summary
  calibration_summary <- NULL
  calibration_by_action <- NULL
  
  if (!is.null(x$objective_calibration)) {
    
    cal <- x$objective_calibration
    
    calibration_summary <- data.frame(
      scale_cost =
        if (!is.null(cal$scale_cost)) cal$scale_cost else NA_real_,
      
      scale_config =
        if (!is.null(cal$scale_config)) cal$scale_config else NA_real_,
      
      scale_target =
        if (!is.null(cal$scale_target)) cal$scale_target else NA_real_,
      
      n_actions_requested =
        if (!is.null(cal$n_actions_requested)) {
          cal$n_actions_requested
        } else {
          NA_integer_
        },
      
      n_actions_sampled =
        if (!is.null(cal$n_actions_sampled)) {
          cal$n_actions_sampled
        } else {
          NA_integer_
        },
      
      n_moves_requested =
        if (!is.null(cal$n_moves_requested)) {
          cal$n_moves_requested
        } else {
          NA_integer_
        },
      
      n_moves_valid =
        if (!is.null(cal$n_moves_valid)) {
          cal$n_moves_valid
        } else {
          NA_integer_
        },
      
      n_cost_nonzero =
        if (!is.null(cal$n_cost_nonzero)) {
          cal$n_cost_nonzero
        } else {
          NA_integer_
        },
      
      n_config_nonzero =
        if (!is.null(cal$n_config_nonzero)) {
          cal$n_config_nonzero
        } else {
          NA_integer_
        },
      
      n_target_nonzero =
        if (!is.null(cal$n_target_nonzero)) {
          cal$n_target_nonzero
        } else {
          NA_integer_
        },
      
      calibrated =
        if (!is.null(cal$calibrated)) cal$calibrated else NA,
      
      stringsAsFactors = FALSE
    )
    
    if (!is.null(cal$by_action)) {
      
      calibration_by_action <- cal$by_action
      
      if (!is.null(action_names) &&
          "action" %in% names(calibration_by_action)) {
        
        action_ids <- calibration_by_action$action
        
        calibration_by_action$action_label <- ifelse(
          action_ids + 1L <= length(action_names),
          action_names[action_ids + 1L],
          paste0("Action ", action_ids)
        )
      }
    }
  }
  
  # Objective preference summary
  objective_summary <- NULL
  
  if (!is.null(x$objective_weights)) {
    
    weights <- x$objective_weights
    
    objective_summary <- data.frame(
      alpha =
        if (!is.null(weights$alpha)) weights$alpha else NA_real_,
      
      beta =
        if (!is.null(weights$beta)) weights$beta else NA_real_,
      
      gamma =
        if (!is.null(weights$gamma)) weights$gamma else NA_real_,
      
      stringsAsFactors = FALSE
    )
  }
  
  # Optimisation-run summary
  run_summary <- data.frame(
    best_objective =
      x$H_best,
    
    iterations_run =
      x$iterations_run,
    
    proposals =
      x$proposals,
    
    accepted =
      x$accepted,
    
    overall_acceptance_percent =
      100 * x$overall_acc,
    
    early_stop_iter =
      x$early_stop_iter,
    
    iteration_time_ms_total =
      if (!is.null(x$iter_ms_total)) x$iter_ms_total else NA_real_,
    
    iteration_time_ms_mean =
      if (!is.null(x$iter_ms_total) &&
          !is.null(x$iter_count) &&
          x$iter_count > 0) {
        x$iter_ms_total / x$iter_count
      } else {
        NA_real_
      },
    
    stringsAsFactors = FALSE
  )
  
  # Interpretable headline statistics
  headline <- list(
    n_species = n_species,
    
    n_species_target_met =
      sum(species_summary$target_met),
    
    prop_species_target_met =
      mean(species_summary$target_met),
    
    mean_target_attainment =
      mean(species_summary$total_target_attainment),
    
    median_target_attainment =
      stats::median(species_summary$total_target_attainment),
    
    mean_target_shortfall =
      mean(species_summary$remaining_target_shortfall),
    
    total_creation_contribution =
      sum(species_summary$creation_target_contribution),
    
    total_improvement_contribution =
      sum(species_summary$improvement_target_contribution),
    
    n_selected_cells =
      if (is.null(action_summary)) {
        NA_integer_
      } else {
        sum(action_summary$selected_cells)
      }
  )
  
  out <- list(
    headline = headline,
    run = run_summary,
    objective = objective_summary,
    calibration = calibration_summary,
    calibration_by_action = calibration_by_action,
    actions = action_summary,
    species = species_summary
  )
  
  class(out) <- "ebrel_summary"
  
  if (isTRUE(print)) {
    print(out)
  }
  
  out
}