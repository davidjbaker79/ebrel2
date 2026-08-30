#' Convert a Flat Vector to a SpatRaster Stack
#'
#' Converts a flattened vector (returned by C++ code in column-major order)
#' into a `SpatRaster` object with multiple layers using the `terra` package.
#'
#' @param vec A numeric vector of length `dim1 * dim2 * nh`, where each segment of length `dim1 * dim2`
#'   corresponds to one habitat layer, stored in column-major order.
#' @param nh Integer. Number of habitat layers (i.e., number of raster bands).
#' @param dim1 Integer. Number of columns in the raster (x dimension).
#' @param dim2 Integer. Number of rows in the raster (y dimension).
#' @param extent Optional. A `SpatExtent` object or extent values (passed to `terra::ext()`).
#'   Defines the spatial extent of the output raster.
#' @param crs Optional. Coordinate reference system string or object (passed to `terra::crs()`).
#'
#' @return A `SpatRaster` object with `nh` layers and dimensions `[dim2, dim1]`.
#'
#' @examples
#' # Example with 2 layers of 3x3 raster
#' vec <- c(1:9, 11:19)
#' r <- vector_to_spatraster(vec, nh = 2, dim1 = 3, dim2 = 3)
#' plot(r)
#'
#' @import terra
#' @export
#'
vector_to_spatraster <- function(vec, nh, dim1, dim2, extent = NULL, crs = NA) {

  # dimensions
  size <- dim1 * dim2
  layers <- vector("list", nh)

  for (h in seq_len(nh)) {
    start <- (h - 1) * size + 1
    end <- h * size

    # C++ stores in column-major order
    layers[[h]]  <- matrix(vec[start:end], nrow = dim2, ncol = dim1, byrow = TRUE)
  }

  # 3D array: [rows, cols, layers] → [dim2, dim1, nh] after transpose
  arr <- array(unlist(layers), dim = c(dim2, dim1, nh))

  # terra expects [rows, cols, layers] → flip dims 1 and 2
  r <- terra::rast(arr)

  if (!is.null(extent)) terra::ext(r) <- extent
  if (!is.na(crs)) terra::crs(r) <- crs

  return(r)
}
