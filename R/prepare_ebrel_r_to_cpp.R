#' Prepare EBREL Inputs for the C++ Backend
#'
#' Converts raster, vector, and matrix inputs into the flattened data structures
#' used by the EBREL C++ backend.
#'
#' The function validates the dimensions and basic coding of the supplied model
#' inputs, extracts raster values, and returns a named list containing flat
#' vectors and model dimensions.
#'
#' @param E_rast A single-layer [terra::SpatRaster] containing the existing
#'   habitat class in each landscape cell. Values must use the internal
#'   zero-indexed habitat coding:
#'   \describe{
#'     \item{`-1`}{No existing habitat.}
#'     \item{`0:(n_habitats - 1)`}{Existing habitat class.}
#'   }
#'
#' @param C_rast A multilayer [terra::SpatRaster] containing one cost layer per
#'   action. The first `n_actions - 1` layers represent habitat creation
#'   actions and the final layer represents habitat improvement. Cell-action
#'   combinations with costs equal or approximately equal to the configured
#'   sentinel value are interpreted as unavailable by downstream C++ code.
#'
#' @param SD_rast A multilayer [terra::SpatRaster] containing one baseline
#'   species-distribution layer per species. Its geometry must match `E_rast`
#'   and `C_rast`.
#'
#' @param D_vec A numeric vector containing one maximum dispersal-distance value
#'   per species. Its length must equal the number of layers in `SD_rast`.
#'
#' @param SxH_mat A numeric, integer, or logical matrix describing
#'   species-habitat associations. Rows correspond to species and columns
#'   correspond to habitat creation classes. It must have dimensions
#'   `n_species` by `n_habitats`, where
#'   `n_habitats = n_actions - 1`.
#'
#' @param O_vec A numeric vector containing one biodiversity target per species.
#'   Its length must equal the number of layers in `SD_rast`.
#'
#' @param sigma A positive numeric scalar controlling the strength of distance
#'   decay used when weighting habitat creation proposals. Proposal weights are
#'   based on
#'   \deqn{w = \exp(-\sigma d),}
#'   where (d) is distance from the relevant species distribution.
#'
#' @param LM_rast A single-layer [terra::SpatRaster] containing the binary land
#'   mask. Values must be `0` or `1`, and its geometry must match the other
#'   raster inputs.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{`dim_x`}{Number of raster columns.}
#'     \item{`dim_y`}{Number of raster rows.}
#'     \item{`n_actions`}{Number of action layers in `C_rast`.}
#'     \item{`n_species`}{Number of species layers in `SD_rast`.}
#'     \item{`n_habitats`}{Number of habitat creation classes, equal to
#'       `n_actions - 1`.}
#'     \item{`E`}{Integer vector of length `dim_x * dim_y`, containing the
#'       categorical existing-habitat map.}
#'     \item{`C`}{Numeric vector of length
#'       `dim_x * dim_y * n_actions`, stored in action-major order.}
#'     \item{`SD`}{Numeric vector of length
#'       `dim_x * dim_y * n_species`, stored in species-major order.}
#'     \item{`D`}{Numeric vector of length `n_species`.}
#'     \item{`SxH`}{Numeric vector of length
#'       `n_species * n_habitats`, stored in species-major row order.}
#'     \item{`O`}{Numeric vector of length `n_species`.}
#'     \item{`sigma`}{Validated distance-decay parameter.}
#'     \item{`LM`}{Integer vector of length `dim_x * dim_y`.}
#'   }
#'
#' @details
#' Raster values from `C_rast` and `SD_rast` are extracted as matrices with
#' cells in rows and layers in columns. Calling `as.vector()` on these matrices
#' produces layer-major vectors: all cells from layer 1 are followed by all
#' cells from layer 2, and so forth. This matches the action-major and
#' species-major layouts used by the C++ backend.
#'
#' `E_rast` and `LM_rast` are single-layer rasters and are therefore returned
#' as vectors of length `dim_x * dim_y`.
#'
#' Because R matrices are stored column-major, `SxH_mat` must be transposed
#' before flattening when the C++ code expects species-major row order:
#' `as.numeric(t(SxH_mat))`.
#'
#' @examples
#' \dontrun{
#' obj <- prepare_ebrel_r_to_cpp(
#'   E_rast = ebrel_sim_data$E,
#'   C_rast = ebrel_sim_data$C,
#'   SD_rast = ebrel_sim_data$SD,
#'   D_vec = ebrel_sim_data$D,
#'   SxH_mat = ebrel_sim_data$SxH,
#'   O_vec = ebrel_sim_data$O,
#'   sigma = ebrel_sim_data$sigma,
#'   LM_rast = ebrel_sim_data$LM
#' )
#' }
#'
#' @export
prepare_ebrel_r_to_cpp <- function(E_rast, 
                                   C_rast, 
                                   SD_rast, 
                                   D_vec, 
                                   SxH_mat, 
                                   O_vec, 
                                   sigma, 
                                   LM_rast,
                                   sentinel = 1e10) {
  
  # ---- Type checks ----
  if (!inherits(E_rast, "SpatRaster")) {
    stop("E_rast must be a terra::SpatRaster.")
  }
  if (!inherits(C_rast, "SpatRaster")) {
    stop("C_rast must be a terra::SpatRaster.")
  }
  if (!inherits(SD_rast, "SpatRaster")) {
    stop("SD_rast must be a terra::SpatRaster.")
  }
  if (!inherits(LM_rast, "SpatRaster")) {
    stop("LM_rast must be a terra::SpatRaster.")
  }
  
  if (!terra::compareGeom(
    E_rast,
    C_rast,
    stopOnError = FALSE
  )) {
    stop("`E_rast` and `C_rast` must have matching geometry.")
  }
  
  if (!terra::compareGeom(
    E_rast,
    SD_rast,
    stopOnError = FALSE
  )) {
    stop("`E_rast` and `SD_rast` must have matching geometry.")
  }
  
  if (!terra::compareGeom(
    E_rast,
    LM_rast,
    stopOnError = FALSE
  )) {
    stop("`E_rast` and `LM_rast` must have matching geometry.")
  }
  
  # ---- Dimensions from rasters ----
  dim_x <- terra::ncol(E_rast)  # x/longitude/easting
  dim_y <- terra::nrow(E_rast)  # y/latitude/northing
  
  n_cells <- terra::ncell(E_rast)
  
  n_actions  <- terra::nlyr(C_rast)
  n_habitats <- n_actions - 1L
  n_species  <- terra::nlyr(SD_rast)
  
  if (terra::nlyr(E_rast) != 1L) {
    stop("E_rast must be a single-layer categorical habitat raster.")
  }
  
  if (terra::nlyr(LM_rast) != 1L) {
    stop("`LM_rast` must contain exactly one land-mask layer.")
  }
  
  if (n_actions < 2L) {
    stop(
      "`C_rast` must contain at least one habitat-creation layer ",
      "and one habitat-improvement layer."
    )
  }
  
  if (n_species < 1L) {
    stop("`SD_rast` must contain at least one species layer.")
  }

  # ---- Scalar parameters ----

  if (
    !is.numeric(sigma) ||
    length(sigma) != 1L ||
    is.na(sigma) ||
    !is.finite(sigma) ||
    sigma <= 0
  ) {
    stop("`sigma` must be a single positive finite numeric value.")
  }
  
  if (
    !is.numeric(sentinel) ||
    length(sentinel) != 1L ||
    is.na(sentinel) ||
    !is.finite(sentinel) ||
    sentinel <= 0
  ) {
    stop("`sentinel` must be a single positive finite numeric value.")
  }
  
  # ---- Species vectors ----
  
  if (
    !is.numeric(D_vec) ||
    is.matrix(D_vec) ||
    is.list(D_vec)
  ) {
    stop("`D_vec` must be a numeric vector.")
  }
  
  if (length(D_vec) != n_species) {
    stop(
      sprintf(
        paste0(
          "`D_vec` length (%d) must equal the number of species ",
          "layers in `SD_rast` (%d)."
        ),
        length(D_vec),
        n_species
      )
    )
  }
  
  if (anyNA(D_vec) || any(!is.finite(D_vec))) {
    stop("`D_vec` must not contain missing or non-finite values.")
  }
  
  if (any(D_vec < 0)) {
    stop("`D_vec` must not contain negative dispersal distances.")
  }
  
  if (
    !is.numeric(O_vec) ||
    is.matrix(O_vec) ||
    is.list(O_vec)
  ) {
    stop("`O_vec` must be a numeric vector.")
  }
  
  if (length(O_vec) != n_species) {
    stop(
      sprintf(
        paste0(
          "`O_vec` length (%d) must equal the number of species ",
          "layers in `SD_rast` (%d)."
        ),
        length(O_vec),
        n_species
      )
    )
  }
  
  if (anyNA(O_vec) || any(!is.finite(O_vec))) {
    stop("`O_vec` must not contain missing or non-finite values.")
  }
  
  if (any(O_vec < 0)) {
    stop("`O_vec` must not contain negative targets.")
  }
  
  # ---- Species-habitat matrix ----
  
  if (!is.matrix(SxH_mat)) {
    stop("`SxH_mat` must be a matrix.")
  }
  
  if (
    nrow(SxH_mat) != n_species ||
    ncol(SxH_mat) != n_habitats
  ) {
    stop(
      sprintf(
        paste0(
          "`SxH_mat` must have dimensions n_species x n_habitats ",
          "(%d x %d); found %d x %d."
        ),
        n_species,
        n_habitats,
        nrow(SxH_mat),
        ncol(SxH_mat)
      )
    )
  }
  
  if (!(
    is.numeric(SxH_mat) ||
    is.integer(SxH_mat) ||
    is.logical(SxH_mat)
  )) {
    stop("`SxH_mat` must be numeric, integer, or logical.")
  }
  
  if (anyNA(SxH_mat)) {
    stop("`SxH_mat` must not contain missing values.")
  }
  
  if (!all(SxH_mat %in% c(0, 1))) {
    stop("`SxH_mat` must contain only 0 and 1.")
  }
  

  # ---- Flatten rasters (layer-major) ----
  E_mat  <- terra::values(E_rast, mat = FALSE)
  C_mat  <- terra::values(C_rast, mat = TRUE)
  SD_mat <- terra::values(SD_rast, mat = TRUE)
  LM_mat <- terra::values(LM_rast, mat = TRUE)

  if (length(E_values) != n_cells) {
    stop("Internal error extracting `E_rast`: unexpected length.")
  }
  
  if (
    !is.matrix(C_matrix) ||
    nrow(C_matrix) != n_cells ||
    ncol(C_matrix) != n_actions
  ) {
    stop("Internal error extracting `C_rast`: unexpected dimensions.")
  }
  
  if (
    !is.matrix(SD_matrix) ||
    nrow(SD_matrix) != n_cells ||
    ncol(SD_matrix) != n_species
  ) {
    stop("Internal error extracting `SD_rast`: unexpected dimensions.")
  }
  
  if (length(LM_values) != n_cells) {
    stop("Internal error extracting `LM_rast`: unexpected length.")
  }

  # ---- Validate E ----
  
  if (anyNA(E_values)) {
    stop("`E_rast` must not contain missing values.")
  }
  
  if (any(E_values != as.integer(E_values))) {
    stop("`E_rast` must contain integer habitat codes.")
  }
  
  E_values <- as.integer(E_values)
  
  if (
    any(
      E_values < -1L |
      E_values >= n_habitats
    )
  ) {
    stop(
      "`E_rast` values must be in the range -1:",
      n_habitats - 1L,
      "."
    )
  }
  
  # ---- Validate C ----
  
  if (anyNA(C_matrix)) {
    stop("`C_rast` must not contain missing values.")
  }
  
  if (any(!is.finite(C_matrix))) {
    stop("`C_rast` must not contain non-finite values.")
  }
  
  if (any(C_matrix < 0)) {
    stop("`C_rast` must not contain negative costs.")
  }
  
  # ---- Validate SD ----
  
  if (anyNA(SD_matrix)) {
    stop("`SD_rast` must not contain missing values.")
  }
  
  if (any(!is.finite(SD_matrix))) {
    stop("`SD_rast` must not contain non-finite values.")
  }
  
  # ---- Validate LM ----
  
  if (anyNA(LM_values)) {
    stop("`LM_rast` must not contain missing values.")
  }
  
  if (any(LM_values != as.integer(LM_values))) {
    stop("`LM_rast` must contain integer values.")
  }
  
  LM_values <- as.integer(LM_values)
  
  if (!all(LM_values %in% c(0L, 1L))) {
    stop("`LM_rast` must contain only 0 and 1.")
  }
  
  # ---- Flatten in C++-aligned order ----
  
  # Action-major:
  # all cells for action 0, then action 1, etc.
  C_values <- as.vector(C_matrix)
  
  # Species-major:
  # all cells for species 0, then species 1, etc.
  SD_values <- as.vector(SD_matrix)
  
  # Species-major row order for the C++ SxH layout.
  SxH_values <- as.numeric(t(SxH_mat))
  
  E  <- as.vector(E_mat)
  C  <- as.vector(C_mat)
  SD <- as.vector(SD_mat)
  LM <- as.vector(LM_mat)

  # ---- Final checks ----
  
  if (length(C_values) != n_cells * n_actions) {
    stop("`C` length mismatch after flattening.")
  }
  
  if (length(SD_values) != n_cells * n_species) {
    stop("`SD` length mismatch after flattening.")
  }
  
  if (length(SxH_values) != n_species * n_habitats) {
    stop("`SxH` length mismatch after flattening.")
  }

  # ---- Return list formatted for C++ use ----
  list(
    dim_x = as.integer(dim_x),
    dim_y = as.integer(dim_y),
    n_actions = as.integer(n_actions),
    n_species = as.integer(n_species),
    n_habitats = as.integer(n_habitats),
    E = E_values,
    C = as.numeric(C_values),
    SD = as.numeric(SD_values),
    D = as.numeric(D_vec),
    SxH = SxH_values,
    O = as.numeric(O_vec),
    sigma = as.numeric(sigma),
    LM = LM_values,
    sentinel = as.numeric(sentinel)
  )
}

