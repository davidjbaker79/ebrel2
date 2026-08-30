#' Plot objective and auxiliary traces at improvement steps
#'
#' Filters a simulated annealing run to the iterations where the objective
#' \eqn{H(x)} achieves a new best (optionally requiring a minimum decrease
#' \code{tol}), and plots the corresponding \code{H} and \code{F} values
#' at those improvement steps.
#'
#' This is useful for diagnosing whether reductions in the main objective
#' \eqn{H} coincide with desirable changes in an auxiliary score \code{F}
#' (e.g., a component of the objective reported by your solver).
#'
#' @param H_trace Numeric vector of objective values \eqn{H(x)} per proposal
#'   (e.g., \code{sa_res$H_trace}).
#' @param F_trace Numeric vector of auxiliary values \code{F(x)} aligned with
#'   \code{H_trace} (same length as \code{H_trace}).
#' @param tol Non-negative numeric scalar. Minimum decrease required for an
#'   iteration to count as a new best (default \code{0} = any strict decrease).
#' @param col_H,col_F Colors for the \code{H} and \code{F} lines.
#' @param legend Logical; add a legend (default \code{TRUE}).
#' @param legend_pos Legend position passed to \code{\link[graphics]{legend}}
#'   (default \code{"topright"}).
#' @param lwd Line width for both series (default \code{2}).
#' @param xlab,ylab Axis labels.
#' @param add Logical; if \code{TRUE}, add to an existing plot instead of
#'   creating a new one.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#'   (when \code{add = FALSE}) or ignored (when \code{add = TRUE}).
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{indices}: integer indices of improvement iterations;
#'   \item \code{Hx}: \code{H_trace[indices]};
#'   \item \code{Fx}: \code{F_trace[indices]}.
#' }
#' If no improvements are found, nothing is plotted and the returned vectors are empty.
#'
#' @examples
#' \dontrun{
#' # Given a run result 'sa_res' with H_trace and F_trace:
#' plot_sa_improvements(sa_res$H_trace, sa_res$F_trace, tol = 1e-6)
#'
#' # Add to an existing plot:
#' plot_sa_improvements(sa_res$H_trace, sa_res$F_trace, col_H = "black",
#'                      col_F = "orange", add = FALSE)
#' }
#' @export
plot_sa_improvements <- function(H_trace, F_trace,
                                 tol = 0,
                                 col_H = "blue",
                                 col_F = "red",
                                 legend = TRUE,
                                 legend_pos = "topright",
                                 lwd = 2,
                                 xlab = "Improvement #",
                                 ylab = "Value",
                                 add = FALSE,
                                 ...) {
  # --- checks ---
  if (!is.numeric(H_trace) || !is.numeric(F_trace))
    stop("`H_trace` and `F_trace` must both be numeric vectors.")
  if (length(H_trace) != length(F_trace))
    stop("`H_trace` and `F_trace` must have the same length.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a single non-negative numeric scalar.")

  n <- length(H_trace)
  if (n == 0L) {
    message("No data: `H_trace` is empty.")
    return(invisible(list(indices = integer(0), Hx = numeric(0), Fx = numeric(0))))
  }

  # --- find improvement indices: H strictly below previous best (by tol) ---
  prior_best <- c(Inf, cummin(H_trace)[-n])  # best up to i-1
  i_imp <- which(H_trace < (prior_best - tol))

  if (length(i_imp) == 0L) {
    message("No improvements found (given tol = ", tol, ").")
    return(invisible(list(indices = integer(0), Hx = numeric(0), Fx = numeric(0))))
  }

  Hx <- H_trace[i_imp]
  Fx <- F_trace[i_imp]
  ix <- seq_along(Hx)

  # --- plot ---
  if (!add) {
    ylim <- range(c(Hx, Fx), na.rm = TRUE)
    plot(ix, Hx, type = "l", lwd = lwd, col = col_H,
         xlab = xlab, ylab = ylab, ylim = ylim, ...)
  } else {
    lines(ix, Hx, lwd = lwd, col = col_H, ...)
  }
  lines(ix, Fx, lwd = lwd, col = col_F, ...)

  if (isTRUE(legend)) {
    graphics::legend(legend_pos,
                     legend = c("H (best-so-far)", "F at improvements"),
                     col = c(col_H, col_F), lty = 1, lwd = lwd, bty = "n")
  }

  invisible(list(indices = i_imp, Hx = Hx, Fx = Fx))
}
