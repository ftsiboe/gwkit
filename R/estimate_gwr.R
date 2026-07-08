# ============================================================
# Geographically weighted local mean/variance regression
# ============================================================
# estimate_gwr_by_point()   - units carry point coordinates.
# estimate_gwr_by_polygon() - units are polygons; centroids from geometry.
#
# For every spatial unit, a geographically weighted regression of
# `response ~ covariate` gives the local slope on the MEAN of the response; a
# second GWR of log(residual^2 + eps) ~ covariate gives the local slope on the
# VARIANCE. Distances and the fixed bandwidth come from GWmodel, but the local
# weighted regression and its slope standard error are computed directly, because
# GWmodel::gwr.basic() returns no standard errors at supplied regression points:
# slope, SE via Var(b) = sigma^2 (X'WX)^-1 (X'W^2X) (X'WX)^-1, t, and a normal-
# approximation p-value. The output schema (mean_/var_ estimate, standard_error,
# t_value, p_value) is agnostic to the spatial-unit identifier.
# ============================================================

utils::globalVariables(c("a_hat", "b_hat", "i.a_hat", "i.b_hat", "log_resid_sq",
                         "longitude", "latitude"))


# ------------------------------------------------------------
# Internal: geographic kernel weights for a distance matrix D (n_obs x n_focal).
# `bw` is a fixed distance (adaptive = FALSE) or a neighbour count.
# ------------------------------------------------------------
.gwr_weights <- function(D, bw, kernel, adaptive) {
  if (adaptive) {
    k   <- min(as.integer(bw), nrow(D))
    bwj <- apply(D, 2L, function(col) sort(col, partial = k)[k])
    U   <- sweep(D, 2L, pmax(bwj, .Machine$double.eps), "/")
  } else {
    U <- D / bw
  }
  switch(kernel,
    gaussian    = exp(-0.5 * U^2),
    exponential = exp(-U),
    bisquare    = ifelse(U < 1, (1 - U^2)^2, 0),
    boxcar      = ifelse(U < 1, 1, 0),
    tricube     = ifelse(U < 1, (1 - U^3)^3, 0),
    stop("Unknown kernel: ", kernel))
}


# ------------------------------------------------------------
# Internal: chunked local weighted regression of y on x at reg points.
# Returns per reg-point estimate (slope), standard_error, t_value, p_value,
# and intercept.
# ------------------------------------------------------------
.gwr_fit <- function(obs_xy, reg_xy, x, y, bw, kernel, adaptive,
                     p, theta, longlat, chunk_size, verbose, label = "") {
  nreg <- nrow(reg_xy)
  est <- se <- tv <- inter <- rep(NA_real_, nreg)
  chunks <- split(seq_len(nreg), ceiling(seq_len(nreg) / chunk_size))
  for (ch in chunks) {
    D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = reg_xy[ch, , drop = FALSE],
                          p = p, theta = theta, longlat = longlat)
    W  <- .gwr_weights(D, bw = bw, kernel = kernel, adaptive = adaptive)
    W2 <- W * W

    Sw    <- colSums(W);          Swx  <- colSums(W * x)
    Swy   <- colSums(W * y);      Swxx <- colSums(W * x^2)
    Swxy  <- colSums(W * x * y);  Swyy <- colSums(W * y^2)
    Sw2   <- colSums(W2);         Sw2x <- colSums(W2 * x)
    Sw2xx <- colSums(W2 * x^2)

    det   <- Sw * Swxx - Swx^2
    slope <- (Sw * Swxy - Swx * Swy) / det
    a     <- (Swy - slope * Swx) / Sw

    # non-negative in exact arithmetic; clamp so floating-point cancellation
    # cannot yield sqrt(negative) = NaN and silently drop near-flat fits.
    SSw    <- pmax(Swyy - 2 * a * Swy - 2 * slope * Swxy +
                   a^2 * Sw + 2 * a * slope * Swx + slope^2 * Swxx, 0)
    sigma2 <- SSw / (Sw - 2)
    var_u  <- pmax(Swx^2 * Sw2 - 2 * Swx * Sw * Sw2x + Sw^2 * Sw2xx, 0) / det^2
    se_b   <- sqrt(sigma2 * var_u)
    tvj    <- slope / se_b

    bad <- !is.finite(det) | det == 0 | (Sw - 2) <= 0 |
           !is.finite(se_b) | se_b <= 0
    slope[bad] <- NA_real_; se_b[bad] <- NA_real_
    tvj[bad]   <- NA_real_; a[bad]    <- NA_real_

    est[ch] <- slope; se[ch] <- se_b; tv[ch] <- tvj; inter[ch] <- a
    if (verbose) message("  gwr ", label, " chunk ", ch[1L], "-",
                         ch[length(ch)], " / ", nreg)
  }
  list(estimate = est, standard_error = se, t_value = tv,
       p_value = 2 * stats::pnorm(-abs(tv)), intercept = inter)
}


# ------------------------------------------------------------
# Internal engine shared by the point and polygon variants. `obs` must already
# carry the coordinate columns named in `coords`.
# ------------------------------------------------------------
.estimate_gwr_core <- function(obs, unit, response, covariate, coords,
                               kernel, adaptive, distance_metric,
                               bw, bw_approach, bw_sample, chunk_size,
                               eps, variance, verbose) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwr_*().")
  if (!requireNamespace("sp", quietly = TRUE))
    stop("Package 'sp' is required for estimate_gwr_*().")

  obs <- data.table::as.data.table(obs)
  dm  <- resolve_distance_metric(distance_metric)
  p <- dm$p; theta <- dm$theta; longlat <- dm$longlat

  keep <- is.finite(obs[[response]]) & is.finite(obs[[covariate]]) &
          is.finite(obs[[coords[1L]]]) & is.finite(obs[[coords[2L]]])
  obs <- obs[keep]
  if (nrow(obs) < 10L) return(NULL)

  obs_xy <- as.matrix(obs[, ..coords])
  reg    <- unique(obs[, c(unit, coords), with = FALSE], by = unit)
  reg_xy <- as.matrix(reg[, ..coords])
  x <- as.numeric(obs[[covariate]])

  # bandwidth selection on a random subsample
  if (is.null(bw)) {
    idx <- sample.int(nrow(obs), min(bw_sample, nrow(obs)))
    sub <- sp::SpatialPointsDataFrame(
      coords = obs_xy[idx, , drop = FALSE],
      data   = data.frame(y = as.numeric(obs[[response]])[idx], x = x[idx]))
    bw <- GWmodel::bw.gwr(y ~ x, data = sub, approach = bw_approach,
                          kernel = kernel, adaptive = adaptive,
                          p = p, theta = theta, longlat = longlat)
  }

  # (1) local MEAN slope
  mean_fit <- .gwr_fit(obs_xy, reg_xy, x = x, y = as.numeric(obs[[response]]),
                       bw = bw, kernel = kernel, adaptive = adaptive,
                       p = p, theta = theta, longlat = longlat,
                       chunk_size = chunk_size, verbose = verbose, label = "mean")

  out <- data.table::data.table(
    term                = covariate,
    unit_id             = as.character(reg[[unit]]),
    unit_level          = unit,
    mean_estimate       = mean_fit$estimate,
    mean_standard_error = mean_fit$standard_error,
    mean_t_value        = mean_fit$t_value,
    mean_p_value        = mean_fit$p_value)

  # (2) local VARIANCE slope, from residuals of each unit's own local mean fit
  if (isTRUE(variance)) {
    reg[, `:=`(a_hat = mean_fit$intercept, b_hat = mean_fit$estimate)]
    obs[reg, on = unit, `:=`(a_hat = i.a_hat, b_hat = i.b_hat)]
    obs[, log_resid_sq :=
          log((as.numeric(get(response)) - (a_hat + b_hat * as.numeric(get(covariate))))^2 + eps)]
    var_fit <- .gwr_fit(obs_xy, reg_xy, x = x, y = obs$log_resid_sq,
                        bw = bw, kernel = kernel, adaptive = adaptive,
                        p = p, theta = theta, longlat = longlat,
                        chunk_size = chunk_size, verbose = verbose, label = "var")
    out[, `:=`(var_estimate       = var_fit$estimate,
               var_standard_error = var_fit$standard_error,
               var_t_value        = var_fit$t_value,
               var_p_value        = var_fit$p_value)]
  }

  out[, model_estimator := "gwr"]
  out[]
}


#' Geographically weighted local mean/variance regression at points
#'
#' For each point-referenced unit, fits a geographically weighted regression of
#' `response` on `covariate` (the local slope on the response mean) and, unless
#' disabled, a second GWR of the log-squared residuals (the local slope on the
#' response variance). Local slope standard errors, t-values, and normal-
#' approximation p-values are computed directly (GWmodel returns no SE at
#' supplied regression points). Distances and the fixed bandwidth use GWmodel.
#'
#' @param data A `data.table`/data frame with one row per observation, containing
#'   the unit id (`unit`), `response`, `covariate`, and the two coordinate
#'   columns named in `coords`. Multiple observations per unit (e.g. a time
#'   series) are expected for a meaningful slope.
#' @param unit Character string naming the spatial-unit id column.
#' @param response Character string naming the response column.
#' @param covariate Character string naming the covariate (e.g. `"trend"`).
#'   Default `"trend"`.
#' @param coords Length-2 character vector naming the longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param kernel GW kernel for `GWmodel::bw.gwr()` and the local weights. One of
#'   `"gaussian"`, `"exponential"`, `"bisquare"`, `"boxcar"`, `"tricube"`.
#'   Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`, fixed distance if
#'   `FALSE`. Default `FALSE`.
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Great Circle"`.
#' @param bw Optional pre-computed bandwidth. If `NULL` it is selected once via
#'   `GWmodel::bw.gwr()` on a random subsample.
#' @param bw_approach Bandwidth criterion for `bw.gwr()`. Default `"AICc"`.
#' @param bw_sample Observations sampled for bandwidth selection. Default `1500`.
#' @param chunk_size Regression points per distance-matrix chunk. Default `250`.
#' @param eps Added to squared residuals before the log. Default `1e-8`.
#' @param variance Logical; also estimate the local variance slope. Default
#'   `TRUE`.
#' @param verbose Logical; print chunk progress. Default `FALSE`.
#'
#' @return A `data.table` with one row per unit: `term`, `unit_id`,
#'   `unit_level`, `mean_estimate`, `mean_standard_error`, `mean_t_value`,
#'   `mean_p_value`, (when `variance`) the `var_*` equivalents, and
#'   `model_estimator` (`"gwr"`).
#' @family Geographically weighted regression
#' @export
estimate_gwr_by_point <- function(data, unit, response, covariate = "trend",
                                  coords = c("longitude", "latitude"),
                                  kernel = "bisquare", adaptive = FALSE,
                                  distance_metric = "Great Circle",
                                  bw = NULL, bw_approach = "AICc",
                                  bw_sample = 1500L, chunk_size = 250L,
                                  eps = 1e-8, variance = TRUE, verbose = FALSE) {
  .estimate_gwr_core(data, unit = unit, response = response, covariate = covariate,
                     coords = coords, kernel = kernel, adaptive = adaptive,
                     distance_metric = distance_metric, bw = bw,
                     bw_approach = bw_approach, bw_sample = bw_sample,
                     chunk_size = chunk_size, eps = eps, variance = variance,
                     verbose = verbose)
}


#' Geographically weighted local mean/variance regression on polygons
#'
#' As `estimate_gwr_by_point()`, but the spatial units are polygons: centroid
#' coordinates are derived from a supplied geometry (`terra` `SpatVector` or
#' `sf`) and joined to `data` by the unit id, then the GWR is run on the
#' centroids.
#'
#' @param data A `data.table`/data frame with one row per observation, containing
#'   the unit id (`unit`), `response`, and `covariate`.
#' @param unit Character string naming the spatial-unit id column in `data`.
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit`.
#' @param response Character string naming the response column.
#' @param covariate Character string naming the covariate. Default `"trend"`.
#' @param poly_id Name of the id field in `polygons`. Default: value of `unit`.
#' @param kernel,adaptive,distance_metric,bw,bw_approach,bw_sample,chunk_size,eps,variance,verbose
#'   See `estimate_gwr_by_point()`.
#'
#' @return A `data.table` with one row per unit (see `estimate_gwr_by_point()`).
#' @family Geographically weighted regression
#' @export
estimate_gwr_by_polygon <- function(data, unit, polygons, response,
                                    covariate = "trend", poly_id = unit,
                                    kernel = "bisquare", adaptive = FALSE,
                                    distance_metric = "Great Circle",
                                    bw = NULL, bw_approach = "AICc",
                                    bw_sample = 1500L, chunk_size = 250L,
                                    eps = 1e-8, variance = TRUE, verbose = FALSE) {

  dm <- resolve_distance_metric(distance_metric)

  cds <- if (inherits(polygons, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required for SpatVector polygons.")
    pg  <- if (isTRUE(dm$longlat) && !terra::is.lonlat(polygons))
             terra::project(polygons, "EPSG:4326") else polygons
    cen <- terra::centroids(pg); xy <- terra::crds(cen)
    data.table::data.table(uid = as.character(terra::values(cen)[[poly_id]]),
                           longitude = xy[, 1L], latitude = xy[, 2L])
  } else if (inherits(polygons, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for sf polygons.")
    pg  <- if (isTRUE(dm$longlat)) sf::st_transform(polygons, 4326) else polygons
    cen <- suppressWarnings(sf::st_centroid(sf::st_geometry(pg)))
    xy  <- sf::st_coordinates(cen)
    data.table::data.table(uid = as.character(sf::st_drop_geometry(pg)[[poly_id]]),
                           longitude = xy[, 1L], latitude = xy[, 2L])
  } else {
    stop("`polygons` must be a terra SpatVector or an sf object.")
  }

  d <- data.table::as.data.table(data)
  d[[unit]] <- as.character(d[[unit]])
  d <- merge(d, cds, by.x = unit, by.y = "uid", all.x = TRUE)

  .estimate_gwr_core(d, unit = unit, response = response, covariate = covariate,
                     coords = c("longitude", "latitude"), kernel = kernel,
                     adaptive = adaptive, distance_metric = distance_metric,
                     bw = bw, bw_approach = bw_approach, bw_sample = bw_sample,
                     chunk_size = chunk_size, eps = eps, variance = variance,
                     verbose = verbose)
}