#' Run an EBREL Optimisation from Terra Rasters
#'
#' High-level user-facing wrapper for running an EBREL spatial optimisation
#' from raster-based model inputs.
#'
#' The function validates and aligns the supplied rasters, extracts their values
#' into the flat vector layouts expected by the C++ implementation, runs the
#' optimisation, and reconstructs the best solution as a
#' [terra::SpatRaster].
#'
#' An optional initial solution may be supplied through `X0`. When `X0` is
#' `NULL`, an initial solution is generated internally by the C++ optimiser.
#'
#' @param E A single-layer [terra::SpatRaster] containing the existing habitat
#'   class in each cell. Values must use the internal zero-indexed habitat
#'   coding:
#'   \describe{
#'     \item{`-1`}{No existing habitat.}
#'     \item{`0:(n_habitats - 1)`}{Existing habitat class.}
#'   }
#'
#' @param C A multilayer [terra::SpatRaster] containing the cost of applying
#'   each action in each cell. It must contain one layer per action. The first
#'   `n_actions - 1` layers represent habitat creation actions and the final
#'   layer represents habitat improvement. Costs approximately equal to
#'   `sentinel` are treated as unavailable actions.
#'
#' @param SD A multilayer [terra::SpatRaster] containing one species
#'   distribution layer per species. Raster geometry must match `C`.
#'
#' @param D A numeric or integer vector containing one dispersal-distance value
#'   per species. Its length must equal the number of layers in `SD`.
#'
#' @param SxH A numeric matrix with species in rows and habitat creation classes
#'   in columns. It must have dimensions `n_species` by `n_habitats`, where
#'   `n_habitats = n_actions - 1`. Non-zero entries indicate that a habitat is
#'   suitable for a species.
#'
#' @param O A numeric vector containing one biodiversity target per species.
#'   Its length must equal the number of layers in `SD`.
#'
#' @param action_weight A numeric vector containing one contribution weight per
#'   action. Its length must equal the number of layers in `C`. The final value
#'   gives the contribution weight assigned to habitat improvement.
#'
#' @param LM A single-layer [terra::SpatRaster] containing the binary land mask.
#'   Values must be `0` or `1`, and raster geometry must match `C`.
#'
#' @param X0 An optional single-layer [terra::SpatRaster] containing the initial
#'   action configuration. Raster geometry must match `C`. Values must use the
#'   internal EBREL action coding:
#'   \describe{
#'     \item{`-1`}{No action selected.}
#'     \item{`0:(n_actions - 1)`}{Zero-indexed selected action.}
#'   }
#'   When `NULL`, an initial solution is generated internally according to the
#'   options supplied in `opt`. Defaults to `NULL`.
#'
#' @param opt A named list of optimisation, calibration, temperature-tuning,
#'   stopping-rule, random-seed, and output options passed to the C++ runner.
#'   Unspecified options use their internal defaults.
#'
#' @param sentinel Numeric scalar identifying unavailable cell-action
#'   combinations in `C`. Defaults to `1e10`.
#'
#' @return A named list containing optimisation results and diagnostics returned
#'   by the C++ runner. The `X_best` element is reconstructed as a single-layer
#'   [terra::SpatRaster] with the same geometry as `E`. Its values use the
#'   internal EBREL action coding:
#'   \describe{
#'     \item{`-1`}{No action selected.}
#'     \item{`0:(n_actions - 1)`}{Zero-indexed selected action.}
#'   }
#'
#'   Other elements may include the initial solution, objective value,
#'   species-target contributions, objective-component traces, acceptance
#'   diagnostics, iteration counts, objective scaling, and temperature-tuning
#'   diagnostics.
#'
#' @details
#' Raster values are flattened in action-major or species-major order as
#' required by the C++ implementation. For `C`, all cells from the first action
#' layer are followed by all cells from the second action layer, and so on.
#' Species distribution layers in `SD` are flattened in the equivalent
#' species-major order.
#'
#' `SxH` is transposed before conversion so that it is passed to C++ in
#' species-major, row-wise order.
#'
#' When `X0` is `NULL`, the wrapper passes an empty integer vector to C++. This
#' signals that the initial action configuration should be generated internally.
#' A supplied `X0` is passed unchanged after validation; it must therefore
#' already use zero-indexed internal action codes.
#'
#' @seealso [hot_start_greedy()]
#'
#' @examples
#' \dontrun{
#' result <- run_ebrel_R(
#'   E = ebrel_sim_data$E,
#'   C = ebrel_sim_data$C,
#'   SD = ebrel_sim_data$SD,
#'   D = ebrel_sim_data$D,
#'   SxH = ebrel_sim_data$SxH,
#'   O = ebrel_sim_data$O,
#'   action_weight = c(
#'     rep(1, terra::nlyr(ebrel_sim_data$C) - 1L),
#'     0.1
#'   ),
#'   LM = ebrel_sim_data$LM,
#'   X0 = NULL,
#'   opt = list(
#'     n_iterations = 10000,
#'     rng_seed = 123,
#'     verbose = TRUE
#'   )
#' )
#'
#' terra::plot(result$X_best)
#'
#' X0 <- hot_start_greedy(ebrel_sim_data)
#'
#' hot_result <- run_ebrel_R(
#'   E = ebrel_sim_data$E,
#'   C = ebrel_sim_data$C,
#'   SD = ebrel_sim_data$SD,
#'   D = ebrel_sim_data$D,
#'   SxH = ebrel_sim_data$SxH,
#'   O = ebrel_sim_data$O,
#'   action_weight = c(
#'     rep(1, terra::nlyr(ebrel_sim_data$C) - 1L),
#'     0.1
#'   ),
#'   LM = ebrel_sim_data$LM,
#'   X0 = X0
#' )
#' }
#'
#' @export
run_ebrel_R <- function(E,
                        C,
                        SD,
                        D,
                        SxH,
                        O,
                        action_weight,
                        LM,
                        X0 = NULL,
                        opt = list(),
                        sentinel = 1e10) {

  # ---- Check alignment of rasters ----
  stopifnot(terra::compareGeom(C, LM, stopOnError = FALSE))
  stopifnot(terra::compareGeom(C, E,  stopOnError = FALSE))
  stopifnot(terra::compareGeom(C, SD, stopOnError = FALSE))

  # ---- Checks ----
  stopifnot(
    inherits(E, "SpatRaster"),
    inherits(C, "SpatRaster"),
    inherits(SD, "SpatRaster"),
    inherits(LM, "SpatRaster")
    )

  # ---- Check X0 ----
  if (!is.null(X0) && !inherits(X0, "SpatRaster")) {
    stop("X0 must be a SpatRaster or NULL")
  }
  
  # ---- Dimensions ----
  dim_x <- terra::ncol(C)
  dim_y <- terra::nrow(C)
  n_cells <- terra::ncell(C)
  n_actions <- terra::nlyr(C)
  n_species <- terra::nlyr(SD)
  n_habitats <- n_actions - 1
    
  # ---- Checks ----
  if (terra::nlyr(E) != 1) {
    stop("E must be a single-layer raster with values -1..(n_habitats)")
  }
  if (length(D) != n_species) {
    stop("D length must equal n_species")
  }
  if (terra::nlyr(LM) != 1L) {
    stop("LM must be a single-layer raster")
  }
  if (!is.null(X0) && terra::nlyr(X0) != 1L) {
    stop("X0 must be a single-layer raster")
  }
  if (!is.matrix(SxH) || nrow(SxH) != n_species || ncol(SxH) != n_habitats) {
    stop("SxH must be n_species x n_habitats")
  }
  if (length(action_weight) != n_actions) {
    stop("action_weight length must equal n_actions")
  }
  if (length(O) != n_species) {
    stop("O length must equal n_species")
  }

  # ---- Template for reconstructed output ----
  tp <- E
  
  # ---- Extract & flatten ---
  E_values  <- as.integer(terra::values(E, mat = FALSE))
  C_values  <- as.vector(terra::values(C, mat = TRUE))
  SD_values <- as.vector(terra::values(SD, mat = TRUE))
  LM_values <- as.integer(terra::values(LM, mat = FALSE))

  # ---- Validate E ----
  if (anyNA(E_values)) {
    stop("E contains NA values")
  }
  
  if (any(E_values < -1L | E_values >= n_habitats)) {
    stop(
      "E values must be in -1:",
      n_habitats - 1L
    )
  }
  
  # ---- Validate LM ----
  if (anyNA(LM_values)) {
    stop("LM contains NA values")
  }
  
  if (!all(LM_values %in% c(0L, 1L))) {
    stop("LM must contain only 0 and 1")
  }

  # ---- Convert optional X0 ----
  if (is.null(X0)) {
    
    # Empty vector signals to C++ that X0 should be generated.
    X0_values <- integer(0)
    
  } else {
    
    X0_values <- terra::values(
      X0,
      mat = FALSE
    )
    
    if (length(X0_values) != n_cells) {
      stop("X0 must contain one value per landscape cell")
    }
    
    if (anyNA(X0_values)) {
      stop("X0 contains NA values")
    }
    
    if (any(X0_values != as.integer(X0_values))) {
      stop("X0 must contain integer action codes")
    }
    
    X0_values <- as.integer(X0_values)
    
    # Internal coding:
    #   -1 = no action
    #    0:(n_actions - 1) = selected action
    if (any(X0_values < -1L | X0_values >= n_actions)) {
      stop(
        "X0 values must be in -1:",
        n_actions - 1L
      )
    }
  }
  
  
  # ---- Call C++ one-shot runner ----
  res <- run_ebrel_cpp(
    E_int  = E_values,
    C  = C_values,
    SD = SD_values,
    D  = as.integer(D),
    SxH = as.numeric(t(SxH)),
    O   = as.numeric(O),
    action_weight = as.numeric(action_weight),
    LM_int  = LM_values,
    X0 = X0_values,
    dim_x = dim_x,
    dim_y = dim_y,
    n_actions   = n_actions,
    n_species = n_species,
    sentinel = sentinel,
    opt = opt
  )

  # --- Reconstruct raster ---
  X_best <- terra::rast(tp)
  terra::values(X_best) <- res$X_best
  names(X_best) <- "X_best"

  res$X_best <- X_best
  res
}
