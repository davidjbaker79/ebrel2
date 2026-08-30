#' @keywords internal
#' @useDynLib ebrel2, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

#' Generate naive starting point
#'
#' @export
generate_X0_CI_R <- function(U, n_actions, dim_x, dim_y, base_prob = 0.85, seed) {
  .Call("_ebrel_generate_X0_CI_R", PACKAGE = "ebrel", U, n_actions, dim_x, dim_y, base_prob, seed)
}

