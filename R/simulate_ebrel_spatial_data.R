#' Simulate spatial habitat, cost, and species distribution data for EBREL
#'
#' Generates a synthetic landscape for testing EBREL workflows, including
#' habitat configuration, land mask, action costs, action availability, species
#' distributions, habitat associations, dispersal distances, and occupancy
#' targets.
#'
#' Habitat layers are first generated as spatially structured rasters using
#' \code{.generate_habitat_rast()}. Species are then assigned to one or more
#' habitat types and seeded into suitable cells. Their distributions are grown
#' outward using a breadth-first search via \code{.bfs_fill()}, constrained by
#' the selected habitat cells and each species' dispersal distance. Dispersal
#' distances are drawn from a Gamma distribution parameterized from
#' \code{disp_longtail} and the maximum grid dimension using
#' \code{.map_longtail_to_gamma()}.
#'
#' The function returns rasters and vectors already arranged for EBREL-style
#' testing in R. Habitat identity is returned as a single integer-coded raster
#' (\code{E}), with \code{-1} indicating no habitat. Cost rasters (\code{C})
#' include one layer per habitat action plus an additional improvement layer.
#' Availability rasters (\code{U}) are derived from the cost rasters, with
#' unavailable actions indicated by 1.
#'
#' @param dim_x Integer. Number of grid rows.
#' @param dim_y Integer. Number of grid columns.
#' @param n_habitats Integer. Number of habitat types. The final habitat type is
#'   treated as unavailable/urban in the cost construction.
#' @param n_species Integer. Number of species to simulate.
#' @param disp_max Integer. Maximum allowed dispersal distance in cells.
#'   Simulated dispersal distances are truncated to this value.
#' @param disp_longtail Integer from 0 to 3. Controls the mean and spread of
#'   species dispersal distances through \code{.map_longtail_to_gamma()}:
#'   \itemize{
#'     \item \code{0}: very short-tailed dispersal
#'     \item \code{1}: short-tailed dispersal
#'     \item \code{2}: intermediate dispersal
#'     \item \code{3}: longer-tailed dispersal
#'   }
#' @param rarity_bias Numeric >= 0. Controls the skew in species prevalence.
#'   Higher values produce more rare species on average.
#' @param fixed_O Optional numeric scalar. If supplied, the same occupancy target
#'   is assigned to all species. If \code{NULL}, species-specific targets are
#'   drawn at random.
#' @param unavail_hab_prop Numeric. Proportion of habitat generation allocated to
#'   unavailable habitat in \code{.generate_habitat_rast()}.
#' @param maskLandArea Logical. If \code{TRUE}, habitat and species rasters are
#'   masked to land cells. If \code{FALSE}, the land mask is still returned but
#'   not applied to the habitat and species layers.
#' @param prop_sea Numeric in \eqn{(0,1)}. Proportion of the grid influenced by
#'   sea when generating the land mask.
#' @param habitat_names Optional character vector of habitat names. If supplied,
#'   these are used to name the cost and availability layers, with
#'   \code{"improvement"} appended as the final action layer.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{dim_x}}{Number of grid rows.}
#'   \item{\code{dim_y}}{Number of grid columns.}
#'   \item{\code{n_actions}}{Number of actions, equal to \code{n_habitats + 1}
#'   because an additional improvement action is included.}
#'   \item{\code{n_species}}{Number of simulated species.}
#'   \item{\code{n_habitats}}{Number of habitat types.}
#'   \item{\code{E}}{A \code{SpatRaster} of integer-coded habitat identity.
#'   Values run from \code{0} to \code{n_habitats - 1}; cells with no habitat
#'   are coded as \code{-1}.}
#'   \item{\code{C}}{A \code{SpatRaster} of action costs with one layer per
#'   habitat action plus one improvement layer. Unavailable actions are assigned
#'   a large sentinel cost.}
#'   \item{\code{U}}{A \code{SpatRaster} indicating action unavailability, derived
#'   from \code{C}. Cells with unavailable actions are coded as 1.}
#'   \item{\code{SD}}{A \code{SpatRaster} with one binary presence-absence layer
#'   per species.}
#'   \item{\code{SxH}}{A binary matrix of dimension
#'   \code{n_species x n_habitats} indicating species-habitat associations.}
#'   \item{\code{D}}{An integer vector of length \code{n_species} giving
#'   species-specific dispersal distances in cells.}
#'   \item{\code{O}}{A numeric vector of length \code{n_species} giving occupancy
#'   targets.}
#'   \item{\code{sigma}}{A scalar equal to \code{1 / mean(D)}.}
#'   \item{\code{LM}}{A binary \code{SpatRaster} land mask.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_ebrel_spatial_data(
#'   dim_x = 50,
#'   dim_y = 50,
#'   n_habitats = 4,
#'   n_species = 15,
#'   disp_max = 20,
#'   rarity_bias = 2,
#'   seed = 123
#' )
#'
#' # Dispersal distances
#' hist(sim$D, main = "Species dispersal distances")
#'
#' # Plot habitat identity raster
#' plot(sim$E, main = "Habitat identity")
#'
#' # Plot first species distribution
#' plot(sim$SD[[1]], main = "Species 1")
#' }
#'
#' @import terra
#' @import gstat
#' @export
simulate_ebrel_spatial_data <- function(
    dim_x = 100,
    dim_y = 100,
    res_m = 100,
    n_habitats = 5, # Number of existing habitats 
    n_species = 5, # Number of species
    disp_max = 25,
    disp_longtail = 0,
    rarity_bias = 1,
    fixed_O = NULL,
    unavail_hab_prop = 0.1,
    maskLandArea = TRUE,
    prop_sea = 0.2,
    habitat_names = NULL,
    seed = NULL) {
  # --- Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # --- Number of cells
  n_cells <- dim_x * dim_y

  # --- Max dispersal
  dim_max <- max(dim_x, dim_y) 

  # --- Generate an autocorrelated response variable
  hab_rast <- .generate_habitat_rast(
    dim_x = dim_x,
    dim_y = dim_y,
    res_m = res_m,
    n_habitats = n_habitats,
    unavail_hab_prop = unavail_hab_prop,
    prop_sea = prop_sea,
    habitat_names = habitat_names
  )

  # --- Mask land or reset landMask
  landMask <- hab_rast$LM
  if (maskLandArea) {
    hab_rast[[2]] <- mask(hab_rast[[2]], landMask)
    hab_rast[[3]] <- mask(hab_rast[[3]], landMask)
  } else {
    # Blank out landMask if no masking occurred
    landMask <- setValues(landMask, NA)
  }

  # --- Generate rarity weights using beta distribution
  rarity_weights <- sort(rbeta(n_species, shape1 = 1, shape2 = 3 * rarity_bias + 0.1), 
                         decreasing = TRUE)

  # --- Simulated species probability distributions (pd)
  pd <- array(0, dim = c(dim_x, dim_y, n_species))

  # ---  Habitat association matrix
  SxH <- matrix(0, nrow = n_species, ncol = n_habitats)
  rownames(SxH) <- paste0("sp_", 1:n_species)
  colnames(SxH) <- paste0("hab_", 1:n_habitats)

  # --- Dispersal paramsdisp_longtail=0
  gamma_params <- .map_longtail_to_gamma(disp_longtail, dim_max)
  long_disp_shape <- gamma_params$shape
  long_disp_scale <- gamma_params$scale
  
  # --- Species dispersal distances ---
  D <- floor(rgamma(n_species, shape = long_disp_shape, scale = long_disp_scale))
  D <- pmax(1, pmin(D, disp_max)) * res_m

  # --- Binary SD array
  SD_rast <- vector("list", n_species)

  # - Try to create species in landscape s =1
  for (s in 1:n_species) {
    npres <- 0
    tries <- 0L

    # cap habitats drawn and align probs
    k_max <- min(4L, n_habitats)

    repeat {
      tries <- tries + 1
      if (tries > 200L) {
        stop(
          "simulate_ebrel_spatial_data(): failed to place species ", s,
          " after 200 attempts (no eligible habitat cells)."
        )
      }

      # draw how many habitats and which ones
      k <- sample.int(k_max, 1)
      hs <- sample.int((n_habitats), k, replace = FALSE)

      # If arable, only consider suitable as near another habitat (edges)
      if (length(hs) == 1 && hs == (n_habitats - 1)) {
        other <- hab_rast$E[[sample(1:(n_habitats-2), 1)]]
        other_edge <- focal(other, w = 3, fun = "max", na.rm = TRUE)
        agr_mod <- hab_rast$E[[(hs)]] * other_edge
        hs_max <- terra::app(agr_mod, "max")
      } else {
        # cells with any of the chosen habitats
        hs_max <- terra::app(hab_rast$E[[hs]], "max")
      }
      
      # sample seed cells ONLY from >0
      ssize <- sample.int(3L, 1)
      cells <- .safe_sample_cells(hs_max, ssize)
      if (length(cells) == 0L) next # try new hs

      # seed raster
      occ_seed_rast <- terra::setValues(hs_max, 0)
      occ_seed_rast[cells] <- 1

      # grow with BFS, constrained by hs_max
      SD_s <- .bfs_fill(hs_max, occ_seed_rast, width = D[s])
      
      # any presences?
      n_pres_cells <- as.numeric(terra::global(SD_s, "sum", na.rm = TRUE)[1, 1])
      if (!is.na(n_pres_cells) && n_pres_cells > 0) {
        SxH[s, hs] <- 1
        SD_rast[[s]] <- SD_s
        npres <- 1
        break
      }
      # else retry with different hs
    }
  }
  SD_rast <- do.call(c, SD_rast)

  # --- Targets
  if (is.null(fixed_O)) {
    O <- round(runif(n_species, min = 0.1, max = 0.3), 1)
  } else {
    O <- rep(fixed_O, n_species)
  }

  # --- Sigma for weights
  sigma <- 1 / mean(D)

  # --- Costs

  # - Set base costs for land management type per pixel
  hab_base_cost <- rev(.approx_doubling(n_habitats, include_zero = FALSE, round_digits = 0))

  # - Set landscape ramp to scale base costs by geographic area
  x_ramp_vec <- seq(0.25, 1, length.out = dim_y)
  geo_C_ramp <- setValues(hab_rast$HO[[1]], rep(x_ramp_vec, dim_x))

  # - Assign costs i = 1
  C_rast <- lapply(seq_len(n_habitats-1), function(i) {
    C_i <- hab_rast$HO[[i]] # Get opportunity
    BC_i <- hab_base_cost[[i]] # Get base cost
    C_i <- ifel(C_i == 1, BC_i, 0)
    C_i <- C_i * geo_C_ramp
    C_i[C_i == 0] <- 1e10
    C_i
  })
  C_rast <- do.call(c, C_rast)

  # Add sentinel cost for urban (i.e. unavailable everywhere)
  C_urb <- setValues(C_rast[[1]], 1e10)
  C_rast <- c(C_rast, C_urb)
  # Name
  names(C_rast) <- paste0("C", 1:(n_habitats))

  # --- Sort out NAs
  E_rast <- ifel(is.na(hab_rast$E), 0, hab_rast$E)
  SD_rast <- ifel(is.na(SD_rast), 0, SD_rast)
  C_rast <- ifel(is.na(C_rast), 1e10, C_rast)
  LM_rast <- ifel(is.na(landMask), 0, landMask)
  
  # --- Apply land masks
  if (maskLandArea) {
    E_rast <- E_rast * LM_rast
    SD_rast <- SD_rast * LM_rast
  }

  # --- E one-hot to id
  idx1 <- app(E_rast, which.max)          # 1..n_habitats
  s    <- app(E_rast, sum, na.rm = TRUE)  # 0..1
  E_id <- idx1 - 1                        # 0..(n_habitats-1)
  E_id <- ifel(s == 0, -1, E_id)          # set no-habitat cells to -1

  # --- Generate improvement eligibility in C
  # Degraded cells: spatially structured subset of all non-urban and arable land 
  # (arable gets improved via AES-type actions)
  # Improvement cost: varies by existing habitat type
  # Improvement allowed only for habitat IDs 1:(n_habitats - 2)
  improvable_ids <- 1:(n_habitats - 2)
  
  I <- classify(
    ifel((E_id + 1) <= (n_habitats - 2), E_id + 1, 0),
    cbind(1:(n_habitats - 2), 0.1 * hab_base_cost[1:(n_habitats - 2)])
  )
  
  I <- ifel(is.na(I) | I == 0, 1e10, I)
  C_rast <- c(C_rast, I)

  # --- U from C
  U_rast <- ifel(C_rast == 1e10, 1, 0)

  # --- Rename is provided
  if (!is.null(habitat_names)){
    names(C_rast) <- c(habitat_names, "improvement")
    names(U_rast) <- c(habitat_names, "improvement")
  }
  
  # --- Return in R formats
  sim_r <- list(
    dim_x = dim_x,
    dim_y = dim_y,
    n_actions = n_habitats + 1,
    n_species = n_species,
    n_habitats = n_habitats,
    E = E_id,
    C = C_rast,
    U = U_rast,
    SD = SD_rast,
    SxH = SxH,
    D = D,
    O = O,
    sigma = sigma,
    LM = LM_rast
  )
  return(sim_r)
}

#' Generate autocorrelated habitat rasters with availability mask
#'
#' Creates a set of \code{SpatRaster} layers representing (i) continuous habitat
#' suitability surfaces (\code{HP}), (ii) a binary one-hot existing habitat stack
#' (\code{E}), and (iii) a binary habitat opportunity stack (\code{HO}) after
#' masking out an "urban/unavailable" fraction of the landscape.
#'
#' Internally, unconditional Gaussian fields are simulated with \pkg{gstat},
#' rescaled to \eqn{[0,1]}, and the lowest \code{unavail_hab_prop} quantile is
#' treated as unavailable. \code{E} is built by one-hot encoding the
#' layer-wise \code{which.max(HP)} and multiplying by a probabilistic
#' conversion mask; the function retries that conversion until every habitat
#' layer in \code{E} has at least one positive cell.
#'
#' @param dim_x Integer. Number of rows in the grid. Default \code{50}.
#' @param dim_y Integer. Number of columns in the grid. Default \code{50}.
#' @param n_habitats Integer. Number of habitat layers to simulate. Default \code{4}.
#' @param unavail_hab_prop Numeric in \eqn{[0,1)}. Proportion of cells treated as
#'   unavailable (e.g., urban) and masked out of \code{HP}/\code{E}/\code{HO}.
#'   Default \code{0.25}.
#'
#' @return A list with three \code{SpatRaster} stacks (each with \code{n_habitats} layers):
#' \describe{
#'   \item{\code{HP}}{Continuous habitat potentials, names \code{hp1..hpN}.}
#'   \item{\code{E}}{Binary existing habitat (one-hot), names \code{E1..EN}.}
#'   \item{\code{HO}}{Binary habitat opportunities after masking, names \code{HO1..HON}.}
#' }
#'
#' @details
#' This function is stochastic. For reproducibility, set R’s seed (e.g., \code{set.seed(1)})
#' before calling. The internal retry loop ensures \emph{every} \code{E} layer has at least
#' one positive cell.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- generate_habitats_rast(dim_x = 30, dim_y = 40, n_habitats = 3, unavail_hab_prop = 0.2)
#' sim$HP; sim$E; sim$HO
#' terra::global(sim$E, "sum", na.rm = TRUE)  # cells per habitat
#' }
#'
#' @importFrom gstat vgm gstat
#' @importFrom terra rast values app setValues which.max ifel global minmax ncell names
#' @noRd
# dim_x = dim_x
# dim_y = dim_y
# n_habitats = 5
# unavail_hab_prop = 0.1
# prop_sea = 0.2
# habitat_names <- c("woodland","heathland","arableMargin","arable","urban")
.generate_habitat_rast <- function(dim_x = 50, dim_y = 50, res_m = 100,
                                   n_habitats = 5, unavail_hab_prop = 0.1, 
                                   prop_sea = 0.2, habitat_names = NULL) {
  
  # --- Create autocorrelated surface
  xy <- expand.grid(x = 1:dim_x, y = 1:dim_y)
  varioMod <- vgm(psill = 1, range = 40, model = 'Exp', nugget = 0.0)
  zDummy <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE,
                  beta = 0, model = varioMod, nmax = 50)

  # -- Simulate surface for "urban"
  urban_sim <- predict(zDummy, newdata = xy, nsim = 1, debug.level = 0)
  urban <- .rescale01(rast(urban_sim, type = "xyz", crs = "EPSG:27700"))
  urban_q <- quantile(values(urban), prob = unavail_hab_prop)
  urban_binary <- ifel(urban < urban_q, 1, NA)

  # -- Simulate surfaces for "habitat"
  habitat_sim <- predict(zDummy, newdata = xy, nsim = (n_habitats-2), debug.level = 0)
  habitat <- rast(habitat_sim, type = "xyz", crs = "EPSG:27700")
  for (x in seq_len(n_habitats-2)) {
    habitat[[x]] <- .rescale01(habitat[[x]])
  }
  names(habitat) <- paste0("hp", 1:(n_habitats-2))
  
  # -- Establish existing habitat
  # - One-hot encode E (i.e. only one habitat per tile)
  habitat_max <- which.max(habitat)
  E_hab <- rast(lapply(seq_len(n_habitats-2), function(k) 1 * (habitat_max == k)))
  E_hab <- mask(E_hab, urban_binary, inverse = TRUE)
  E_hab[is.na(E_hab)] <- 0

  # - Convert HO to existing habitat by converting land
  layers_zero <- 0
  while (layers_zero == 0) {

    # Simulate converted land (e.g. agriculture)
    convert_sim <- predict(zDummy, newdata = xy, nsim = 1, debug.level = 0)
    convert <- (.rescale01(rast(convert_sim, type = "xyz", crs = "EPSG:27700"))^3) * 0.5

    # Convert to binary (presence/absence) and mask out of E
    convert_binary <- .prob_to_bin(convert)
    E <- E_hab * convert_binary

    # Check that all habitats exist
    # per-layer max (vector, one per layer)
    mx <- terra::global(E, "max", na.rm = TRUE)[, 1]
    zero_max_idx <- which(is.na(mx) | mx == 0)
    if (length(zero_max_idx) == 0)
      layers_zero <- 1
  }
  
  # Arable matrix as explicit layer + add in urban
  arable <- 1 - max(E)   # cells not assigned to any semi-natural type
  E <- c(E, arable) 

  # - Mask urban
  E <- mask(E, urban_binary, inverse = TRUE)
  
  # - Add back in urban (no opportunity, but can be existing)
  E <- c(E, urban_binary)
  E[is.na(E)] <- 0
  
  # Rename
  names(E) <- paste0("E", 1:nlyr(E))
  E[is.na(E)] <- 0

  # Generate HO from latent suitability, masked to arable cells
  HO <- rast(lapply(seq_len(n_habitats-2), function(h) {
    ifel(habitat[[h]] >= 0.5, 1, 0)
  }))
  HO <- mask(HO, urban_binary, inverse = TRUE)

  # HO for urban and arable exists (might also be improved_grassland in the UK)
  # must be included for dimensional consistency, but all sentinel
  # Exclusion from actions handled via U_creation and U_improve
  HO_sentinel <- setValues(habitat[[1]], 0)
  HO <- c(HO, HO_sentinel, HO_sentinel)
  names(HO) <- paste0("HO", 1:nlyr(HO))  
  HO[is.na(HO)] <- 0
  
  # - Make sea mask
  LM <- setValues(habitat[[1]], 0)
  LM <- .make_land_mask(LM, prop_sea = prop_sea, smooth_window = 5)

  # - Rename is provided
  if (!is.null(habitat_names)){
    names(E) <- habitat_names
    names(HO) <- habitat_names
  }

  # - Set crs and res
  ext(habitat) <- c(0, dim_x * res_m, 0, dim_y * res_m)
  crs(habitat) <- "EPSG:27700"
  ext(E) <- c(0, dim_x * res_m, 0, dim_y * res_m)
  crs(E) <- "EPSG:27700"
  ext(HO) <- c(0, dim_x * res_m, 0, dim_y * res_m)
  crs(HO) <- "EPSG:27700"
  ext(LM) <- c(0, dim_x * res_m, 0, dim_y * res_m)
  crs(LM) <- "EPSG:27700"
  
  
  list(
    HP = habitat,
    E = E,
    HO = HO,
    LM = LM
  )
}

#' Generate a simple land–sea mask with a smoothed coastline
#'
#' Creates a single-layer \code{SpatRaster} mask distinguishing land from sea,
#' based on a template raster. The coastline is generated by drawing a random
#' shoreline height (in rows from the bottom of the raster) per column,
#' smoothing it with a moving-average filter, and then filling cells below the
#' coastline as sea.
#'
#' The mask is aligned to the geometry (extent, resolution, CRS) of the input
#' raster \code{r}. Cell values are set to \code{1} on land and \code{NA} over
#' sea.
#'
#' @param r A \code{terra::SpatRaster} object providing the target geometry
#'   (extent, resolution, projection, and dimensions) for the mask.
#' @param prop_sea Numeric in \eqn{(0, 1]}. Approximate maximum proportion of
#'   the raster height (from the bottom) that can be occupied by sea. This
#'   determines the upper bound on the randomly drawn coastline heights.
#' @param smooth_window Integer \eqn{\ge 1}. Width of the moving-average filter
#'   used to smooth the randomly drawn coastline heights across columns. Larger
#'   values produce smoother coastlines.
#' @param seed Optional integer. If provided, passed to \code{set.seed()} to
#'   make the coastline generation reproducible.
#'
#' @return A single-layer \code{terra::SpatRaster} with the same geometry as
#'   \code{r}, named \code{"land"}, where:
#'   \itemize{
#'     \item \code{1} indicates land.
#'     \item \code{NA} indicates sea.
#'   }
#'
#' @details
#' Let \eqn{n_r} and \eqn{n_c} be the number of rows and columns in \code{r}.
#' The function:
#' \enumerate{
#'   \item Computes \code{max_sea_rows = round(nr * prop_sea)}, the maximum
#'     number of rows from the bottom that may be sea.
#'   \item Draws a random coastline height (in rows from the bottom) for each
#'     column, uniformly between \code{1} and \code{max_sea_rows}.
#'   \item Smooths these heights across columns using a moving-average filter
#'     of width \code{smooth_window}.
#'   \item For each column, marks all cells below the smoothed coastline as
#'     sea (\code{NA}) and the remaining cells as land (\code{1}).
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Example template raster (e.g. 100 x 100 grid)
#' r <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#'
#' # Generate a land–sea mask with ~30% sea at the bottom and a smooth coastline
#' land_mask <- make_land_mask(r, prop_sea = 0.3, smooth_window = 7, seed = 42)
#'
#' plot(land_mask)
#' }
#'
#' @keywords internal
#' @noRd
.make_land_mask <- function(r, prop_sea = 0.3, smooth_window = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  nr <- nrow(r)
  nc <- ncol(r)

  # Target maximum sea height (in rows from the bottom)
  max_sea_rows <- max(1, round(nr * prop_sea))

  # Random coastline height per column (from 1 to max_sea_rows)
  coast_raw <- sample(1:max_sea_rows, nc, replace = TRUE)

  # Smooth coastline
  k <- rep(1 / smooth_window, smooth_window)
  coast_smooth <- stats::filter(coast_raw, k, sides = 2)

  # replace NAs at edges with original values
  coast_smooth[is.na(coast_smooth)] <- coast_raw[is.na(coast_smooth)]
  coast <- pmax(1, pmin(max_sea_rows, round(coast_smooth)))

  # Build a matrix where 1 = land, NA = sea
  mat <- matrix(1, nrow = nr, ncol = nc)
  for (j in seq_len(nc)) {
    sea_rows <- (nr - coast[j] + 1):nr   # rows that should be NA for column j
    mat[sea_rows, j] <- NA               # row first, then column
  }

  # Create SpatRaster with same geometry as input
  land_r <- r
  values(land_r) <- as.vector(t(mat))
  names(land_r) <- "land"

  return(land_r)
}


#' Buffer binary cells on a raster
#'
#' Internal function that applies a square buffer of radius \code{width} and
#' returns a spatraster where any neighbor within the window is converted to 1.
#'
#' @param r A \code{SpatRaster} (single layer) with binary values 0/1.
#' @param width Integer \eqn{\ge 0}. Neighborhood radius; window size is
#'   \code{(2*width + 1) x (2*width + 1)}. Default \code{2}.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 after function
#'
#' @keywords internal
#' @noRd
#' @importFrom terra focal ifel
.buffer_cells <- function(r, width = 2) {
  w <- matrix(1, nrow = 2*width + 1, ncol = 2*width + 1)
  out <- focal(r, w = w, fun = max, na.policy = "omit", pad = TRUE, padValue = 0)
  out <- ifel(out > 0, 1, 0)
}

#' BFS-style region growth constrained by suitability and hop width
#'
#' Internal function, that, tarting from seed cells (\code{occ_seed_rast > 0}),
#' repeatedly grows the visited set using \code{buffer_cells()} with radius
#' \code{width}, but only admits cells where \code{hs_max > 0}. Stops when no
#' new cells are added.
#'
#' @param hs_max \code{SpatRaster} (single layer) with 0/1 suitability mask.
#' @param occ_seed_rast \code{SpatRaster} (single layer) with 0/1 initial seeds.
#' @param width Integer \eqn{\ge 0}. Maximum hop size per iteration. Default \code{2}.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 indicating the filled region.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra ifel global
.bfs_fill <- function(hs_max, occ_seed_rast, width = 2, erosion_prob = 0.3) {
  max_steps <- sample(1:10, 1)
  step <- 0
  
  hs <- ifel(hs_max > 0, 1, 0)
  visited <- ifel(occ_seed_rast > 0, 1, 0)
  
  # Random directional bias — gives elongated, non-square shapes
  angle <- runif(1, 0, 2 * pi)
  xy_coords <- xyFromCell(hs_max, 1:ncell(hs_max))
  direction_vals <- cos(angle) * xy_coords[, 1] + sin(angle) * xy_coords[, 2]
  direction_vals <- (direction_vals - min(direction_vals, na.rm = TRUE)) /
    (max(direction_vals, na.rm = TRUE) - min(direction_vals, na.rm = TRUE))
  direction_bias <- setValues(hs_max, direction_vals)

  repeat {
    if (step > max_steps) break
    expanded <- .buffer_cells_circular(visited, width = width)
    candidates <- ifel((expanded == 1) & (hs == 1), 1, 0)
    new <- ifel((candidates == 1) & (visited == 0), 1, 0)

    new_vals <- values(new, mat = FALSE)
    accept_vals <- values(direction_bias * hs_max * runif(1, 0.5, 0.95), mat = FALSE)
    
    new_vals[is.na(new_vals)] <- 0
    accept_vals[is.na(accept_vals)] <- 0
    accept_vals <- pmin(pmax(accept_vals, 0), 1)
    
    idx <- which(new_vals == 1)
    if (length(idx) == 0) break
    new_vals[idx] <- rbinom(length(idx), size = 1, prob = accept_vals[idx])
    values(new) <- new_vals
    
    n_new <- sum(new_vals == 1, na.rm = TRUE)
    if (is.na(n_new) || n_new == 0) break
    visited <- ifel(new == 1, 1, visited)
    step <- step + 1
  }
  
  visited <- .erode_boundary(visited, erosion_prob = erosion_prob)
  
  values(visited) <- rbinom(ncell(visited), size = 1, 
                            prob = values(visited) / 2)
  visited

}

.buffer_cells_circular <- function(rast, width = 2) {
  # Circular focal window
  fw <- focalMat(rast, d = width, type = "circle")
  fw[fw > 0] <- 1
  focal(rast, w = fw, fun = "max", na.rm = TRUE)
}

.erode_boundary <- function(visited, erosion_prob = 0.4) {
  # Identify boundary cells — visited cells adjacent to unvisited
  neighbour_sum <- focal(visited, w = matrix(1, 3, 3), fun = "sum", na.rm = TRUE)
  is_boundary <- ifel(visited == 1 & neighbour_sum < 9, 1, 0)
  
  # Randomly erode boundary cells
  boundary_vals <- values(is_boundary, mat = FALSE)
  boundary_vals[is.na(boundary_vals)] <- 0
  
  idx <- which(boundary_vals == 1)
  if (length(idx) > 0) {
    erode <- rbinom(length(idx), size = 1, prob = erosion_prob)
    boundary_vals[idx] <- erode
  }
  
  # Put back into raster and remove eroded cells from visited
  is_boundary_eroded <- setValues(is_boundary, boundary_vals)
  ifel(is_boundary_eroded == 1, 0, visited)
}

#' Rescale to \eqn{[0,1]}
#'
#' Linearly rescales a \code{SpatRaster} to the unit interval using its layer-wise
#' minimum and maximum; safeguards against zero range. Internal function
#'
#' @param x A \code{SpatRaster} (single layer).
#' @return A \code{SpatRaster} on \eqn{[0,1]}.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra minmax
.rescale01 <- function(x) {
  rng <- terra::minmax(x); (x - rng[1]) / max(1e-9, (rng[2] - rng[1]))
}

#' Bernoulli thresholding of a probability raster
#'
#' Converts a probability raster \code{p} into a binary raster by drawing
#' \code{U(0,1)} per cell and returning \code{1*(U < p)}. Internal function.
#'
#' @param prob_rast A \code{SpatRaster} (single layer) with values in \eqn{[0,1]}.
#' @param seed Optional integer seed; if provided, sets R’s RNG state before sampling.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 with name \code{"binary"}.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra rast values ncell
.prob_to_bin <- function(prob_rast, seed = NULL) {

  # draw uniform(0,1) for each cell
  rnd <- rast(prob_rast)
  values(rnd) <- runif(ncell(rnd))

  out <- (rnd < prob_rast) * 1
  names(out) <- "binary"
  out
}

#' Generate an approximately doubling cost scale
#'
#' Creates a short geometric sequence that ends at \code{max_val} and (optionally)
#' starts with 0. Useful for setting tiered costs where each successive tier is
#' about twice the previous one.
#'
#' @param n Integer (>= 1). Number of values to return.
#' @param max_val Numeric (> 0). The maximum (final) positive value in the sequence.
#' @param include_zero Logical. If \code{TRUE}, the first element is 0 and the
#'   remaining \code{n-1} values form a halving series up to \code{max_val}.
#'   If \code{FALSE}, all \code{n} values are positive.
#' @param round_digits Integer or \code{NULL}. If not \code{NULL}, round the
#'   resulting values to this many decimal places.
#'
#' @return A numeric vector of length \code{n}.
#' \itemize{
#'   \item If \code{include_zero = TRUE} and \code{n > 1}:
#'         \code{c(0, max_val/2^{m-1}, ..., max_val/2, max_val)} with \code{m = n-1}.
#'   \item If \code{include_zero = TRUE} and \code{n = 1}: \code{0}.
#'   \item If \code{include_zero = FALSE}:
#'         \code{max_val/2^{n-1}}, \code{...}, \code{max_val/2}, \code{max_val}.
#' }
#'
#' @examples
#' approx_doubling(4, max_val = 800, include_zero = TRUE)
#' # 0 200 400 800
#'
#' approx_doubling(4, max_val = 800, include_zero = FALSE)
#' # 100 200 400 800
#'
#' approx_doubling(5, max_val = 1000, include_zero = TRUE, round_digits = 0)
#'
#' @keywords internal
#' @noRd
.approx_doubling <- function(n, max_val = 1000, include_zero = TRUE, round_digits = NULL) {
  stopifnot(n >= 1)

  if (include_zero) {
    if (n == 1) return(0)
    m <- n - 1                 # number of positive values
    vals <- max_val / 2^((m - 1):0)  # geometric, ends at max_val
    out <- c(0, vals)
  } else {
    m <- n
    vals <- max_val / 2^((m - 1):0)
    out <- vals
  }

  if (!is.null(round_digits)) out <- round(out, round_digits)
  out
}

#' Map a long-tail dispersal class to Gamma parameters (grid-scaled)
#'
#' Converts a discrete long-tail class into parameters of a
#' \eqn{\mathrm{Gamma}(\text{shape}, \text{scale})} distribution used to draw
#' dispersal distances. The mapping is defined in terms of the maximum possible
#' distance on the grid (\code{max_dist}), allowing dispersal distances to scale
#' automatically with grid size (e.g., 25×25 vs 200×200).
#'
#' Each \code{longtail} class corresponds to a pair of fractions
#' (\code{mean_frac}, \code{sd_frac}) that determine the mean and standard
#' deviation of the Gamma distribution as:
#' \deqn{\mu = \text{mean\_frac} \times \text{max\_dist}}
#' \deqn{\sigma = \text{sd\_frac} \times \text{max\_dist}}
#'
#' These are then converted to Gamma parameters via:
#' \deqn{\text{shape} = (\mu^2) / \sigma^2}
#' \deqn{\text{scale} = \sigma^2 / \mu}
#'
#' Larger \code{longtail} values correspond to increasingly long-tailed
#' dispersal distributions (larger mean fraction, larger SD fraction).
#'
#' @param longtail Integer (0–3). Dispersal class where:
#'   \itemize{
#'     \item \code{0}: very short-tailed (\code{mean_frac = 0.05}, \code{sd_frac = 0.10})
#'     \item \code{1}: short-tailed    (\code{mean_frac = 0.15}, \code{sd_frac = 0.15})
#'     \item \code{2}: intermediate    (\code{mean_frac = 0.25}, \code{sd_frac = 0.25})
#'     \item \code{3}: long-tailed     (\code{mean_frac = 0.35}, \code{sd_frac = 0.35})
#'   }
#'
#' @param max_dist Numeric. The maximum possible cell-to-cell distance on the grid
#'   (typically the grid diagonal, e.g., \code{sqrt(nx^2 + ny^2)}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{shape}: Gamma shape parameter.
#'   \item \code{scale}: Gamma scale parameter.
#' }
#'
#' @examples
#' # Suppose a 100x100 grid:
#' max_dist <- sqrt(100^2 + 100^2)
#'
#' .map_longtail_to_gamma(0, max_dist)
#' .map_longtail_to_gamma(1, max_dist)
#' .map_longtail_to_gamma(3, max_dist)
#'
#' @keywords internal
#' @noRd
.map_longtail_to_gamma <- function(longtail, max_dist) {
  # choose relative behaviour by 'longtail'
  base <- switch(as.character(longtail),
                 "0" = list(mean_frac = 0.05, sd_frac = 0.1),
                 "1" = list(mean_frac = 0.15, sd_frac = 0.15),
                 "2" = list(mean_frac = 0.25, sd_frac = 0.25),
                 "3" = list(mean_frac = 0.35, sd_frac = 0.35),
                 stop("Unknown longtail value")
  )

  mean_D <- base$mean_frac * max_dist
  sd_D   <- base$sd_frac   * max_dist

  shape <- (mean_D^2) / (sd_D^2)
  scale <- (sd_D^2) / mean_D

  list(shape = shape, scale = scale)
}


#' Sample eligible cell indices from a mask raster
#'
#' Returns uniformly sampled cell indices from the set of cells where the mask
#' has finite values greater than 0. If there are fewer eligible cells than
#' requested, the sample size is reduced accordingly. If no cells are eligible,
#' an empty integer vector is returned.
#'
#' @param mask_rast A `terra::SpatRaster` mask; cells with `values(mask_rast) > 0`
#'   (and finite) are considered eligible.
#' @param size Integer (>= 1). Desired number of cells to sample.
#'
#' @return An integer vector of 1-based cell indices (length `<= size`).
#'
#' @examples
#' \dontrun{
#' r <- terra::rast(nrows = 5, ncols = 5)
#' terra::values(r) <- c(0,1,2,NA,0)[(1:25 - 1) %% 5 + 1]
#' safe_sample_cells(r, size = 3)  # indices where values > 0
#' }
#'
#' @importFrom terra values
#' @keywords internal
#' @noRd
.safe_sample_cells <- function(mask_rast, size) {
  v <- terra::values(mask_rast)
  idx <- which(is.finite(v) & v > 0)
  if (length(idx) == 0L) return(integer(0))
  size <- min(size, length(idx))
  sample(idx, size)
}

