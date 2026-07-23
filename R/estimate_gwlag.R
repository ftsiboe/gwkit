# ============================================================
# Geographically weighted spatial lag (neighbour-weighted mean)
# ============================================================
# estimate_gwlag() - one entry; point vs polygon detected by class via
#                    .resolve_gw_geometry() (helpers.R).
#
# For every target unit i, the spatial lag of a value column z is the
# geographically weighted mean of z over the source units j:
#
#     lag_i = sum_j w_ij z_j / sum_j w_ij ,
#
# with weights w_ij from GWmodel::gw.weight() over a GWmodel::gw.dist() distance
# matrix (the SAME routines used elsewhere in the pipeline, so results are
# comparable). When `include_self = FALSE` (the default) the focal unit is
# removed from its own neighbourhood (w_ii := 0), giving a strict neighbour mean
# - e.g. the "neighbouring production" m(z_i) term of an availability metric
# z_i + m(z_i). Multiple value columns are lagged in one pass with a shared
# weight matrix.
# ============================================================

# ------------------------------------------------------------
# Internal engine. `obs` carries the source coordinates (.gw_xo, .gw_yo) and
# value_cols; `tgt` carries the target coordinates (.gw_xt, .gw_yt) and ids.
# ------------------------------------------------------------
.estimate_gwlag_core <- function(obs, tgt, unit, value_cols,
                                 p, theta, longlat, kernel, adaptive, bw,
                                 bw_response, bw_approach, include_self) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwlag_*().")

  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)

  obs_xy <- as.matrix(obs[, c(".gw_xo", ".gw_yo")])
  tgt_xy <- as.matrix(tgt[, c(".gw_xt", ".gw_yt")])

  # bandwidth (once) if not supplied: bw.gwr(bw_response ~ 1) on the source set,
  # exactly as the availability step does (fixed-distance CV by default).
  if (is.null(bw)) {
    if (is.null(bw_response)) bw_response <- value_cols[1L]
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(
      coords = obs_xy,
      data   = data.frame(.y = as.numeric(obs[[bw_response]])))
    bw <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                          kernel = kernel, adaptive = adaptive,
                          p = p, theta = theta, longlat = longlat, dMat = d_self)
    rm(d_self)
  }

  # source x target distance and weights
  #
  # .gw_dist_oriented(), NOT gw.dist() directly: GWmodel returns dp x rp on the
  # Minkowski path but rp x dp when longlat = TRUE. The bare call here had the
  # same defect the GWR engine did - gw.weight() applies an adaptive bandwidth
  # column-wise, so a transposed matrix took the neighbourhood along the wrong
  # axis, and the matrix() reshape below then scrambled the source-to-target
  # mapping. Every Great Circle result from this function was invalid; other
  # metrics were unaffected. See .gw_dist_oriented() in gw_local_engine.R.
  D <- .gw_dist_oriented(obs_xy, tgt_xy, p, theta, longlat)
  W <- GWmodel::gw.weight(vdist = D, bw = bw, kernel = kernel, adaptive = adaptive)
  W <- matrix(W, nrow = nrow(obs_xy), ncol = nrow(tgt_xy))   # obs x target

  # strict neighbour mean: zero a source's weight onto its own target
  if (!isTRUE(include_self)) {
    same <- outer(obs[[unit]], tgt[[unit]], FUN = "==")
    W[same] <- 0
  }

  # All value columns in ONE BLAS pass. crossprod(W, V)[t, c] = sum_o W[o,t]*V[o,c],
  # i.e. exactly the per-column colSums(W * z) this replaces, but O(n^2 * ncol) in
  # BLAS rather than an R loop - the difference that makes many-column lags viable.
  W[!is.finite(W)] <- 0
  denom <- colSums(W)
  V <- as.matrix(as.data.frame(obs)[, value_cols, drop = FALSE])
  storage.mode(V) <- "double"; V[!is.finite(V)] <- 0
  num  <- crossprod(W, V)                                  # target x value_col
  bad  <- !is.finite(denom) | denom == 0
  out <- data.table::data.table(uid = tgt[[unit]])
  for (i in seq_along(value_cols)) {
    lag <- num[, i] / denom
    lag[bad] <- NA_real_
    out[[paste0(value_cols[i], "_LM")]] <- lag
  }
  data.table::setnames(out, "uid", unit)
  attr(out, "bandwidth") <- bw
  out[]
}


#' Geographically weighted spatial lag (neighbour-weighted mean)
#'
#' For each target unit, computes the geographically weighted mean of one or more
#' value columns over the source units, using `GWmodel::gw.weight()` weights on a
#' `GWmodel::gw.dist()` distance matrix. With `include_self = FALSE` the focal
#' unit is excluded from its own neighbourhood, yielding a strict neighbour mean
#' (a spatial lag). Point vs polygon mode is class-detected: pass an `sf`/
#' `SpatVector` (as `data`, or via `geometry` for panels) for polygon units,
#' each reduced to its point-on-surface; otherwise the `coords` columns of `data`
#' are used.
#'
#' @param data A `data.table`/data frame of the SOURCE units (one row per unit)
#'   with `unit`, all `value_cols`, and - in point mode - the `coords` columns;
#'   or an `sf`/`SpatVector` carrying both attributes and geometry (polygon mode).
#' @param unit Character; the unit-id column.
#' @param value_cols Character vector of numeric columns to lag. Each returns a
#'   column named `<value_col>_LM`.
#' @param geometry Optional `sf`/`SpatVector` polygon layer (id field `poly_id`
#'   matching `unit`) supplying per-unit geometry. If `NULL` and `data` is not
#'   spatial, point mode uses `coords`.
#' @param coords Length-2 character vector naming longitude/latitude columns
#'   (point mode). Default `c("longitude", "latitude")`.
#' @param poly_id Name of the id field in the polygon layer. Default: `unit`.
#' @param predict Optional target restriction: a data frame of point targets
#'   (point mode) or a character id vector / data frame with `unit` (polygon
#'   mode). Defaults to the units in `data`.
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `FALSE`
#'   (fixed-distance); note `estimate_gwr()`/`estimate_gwss()` default to `TRUE`.
#' @param bw Optional pre-computed bandwidth. If `NULL`, selected once via
#'   `GWmodel::bw.gwr(bw_response ~ 1)` on the source set.
#' @param bw_response Column used for bandwidth selection when `bw` is `NULL`.
#'   Default: the first `value_cols`.
#' @param bw_approach Bandwidth criterion. Default `"CV"`.
#' @param include_self Logical; keep the focal unit in its own neighbourhood.
#'   Default `FALSE` (strict neighbour mean).
#'
#' @return A `data.table` with one row per target unit: the `unit` id and one
#'   `<value_col>_LM` column per value column. The selected bandwidth is attached
#'   as `attr(., "bandwidth")`.
#' @family Geographically weighted summaries
#' @export
estimate_gwlag <- function(data, unit, value_cols, geometry = NULL,
                           coords = c("longitude", "latitude"), poly_id = unit,
                           predict = NULL, distance_metric = "Euclidean",
                           kernel = "bisquare", adaptive = FALSE,
                           bw = NULL, bw_response = NULL, bw_approach = "CV",
                           include_self = FALSE) {
  dm  <- resolve_distance_metric(distance_metric)
  geo <- .resolve_gw_geometry(data, unit = unit, coords = coords,
                              geometry = geometry, poly_id = poly_id,
                              longlat = dm$longlat, predict = predict)
  .estimate_gwlag_core(obs = geo$obs, tgt = geo$tgt, unit = unit,
                       value_cols = value_cols, p = dm$p, theta = dm$theta,
                       longlat = dm$longlat, kernel = kernel, adaptive = adaptive,
                       bw = bw, bw_response = bw_response, bw_approach = bw_approach,
                       include_self = include_self)
}


# ------------------------------------------------------------
# Multi-kernel companion of .estimate_gwlag_core(): builds the (kernel-independent)
# obs->target distance matrix ONCE and reuses it across kernels, and selects the
# bandwidth once (shared, or once per kernel off the SAME self-distance). Stacks
# the per-kernel lags with a `kernel` column. Memory profile matches the single
# call (one obs x target distance held at a time); the saving is the eliminated
# redundant gw.dist()/bw.gwr() over the 5 kernels of a distance metric.
# ------------------------------------------------------------
.estimate_gwlag_kernels_core <- function(obs, tgt, unit, value_cols,
                                         p, theta, longlat, kernels, adaptive,
                                         bw, bandwidth, bw_response, bw_approach, include_self) {
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwlag_kernels().")
  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)
  obs_xy <- as.matrix(obs[, c(".gw_xo", ".gw_yo")])
  tgt_xy <- as.matrix(tgt[, c(".gw_xt", ".gw_yt")])
  nk <- length(kernels)

  # --- bandwidth(s), keyed by kernel, computed once (self-distance built once) ---
  if (is.null(bw)) {
    if (is.null(bw_response)) bw_response <- value_cols[1L]
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(coords = obs_xy,
                                      data = data.frame(.y = as.numeric(obs[[bw_response]])))
    if (identical(bandwidth, "shared")) {
      b0 <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                            kernel = kernels[[1L]], adaptive = adaptive,
                            p = p, theta = theta, longlat = longlat, dMat = d_self)
      bw_named <- stats::setNames(rep(b0, nk), kernels)
    } else {
      bw_named <- stats::setNames(vapply(kernels, function(km)
        GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                        kernel = km, adaptive = adaptive,
                        p = p, theta = theta, longlat = longlat, dMat = d_self),
        numeric(1)), kernels)
    }
    rm(d_self)
  } else if (is.null(names(bw))) {
    bw_named <- stats::setNames(rep(as.numeric(bw), length.out = nk), kernels)
  } else {
    bw_named <- bw[kernels]
    if (any(is.na(bw_named)))
      stop("named `bw` is missing entries for kernel(s): ",
           paste(kernels[is.na(bw_named)], collapse = ", "))
  }

  # --- obs -> target distances built ONCE, reused across kernels ---
  # Orientation-normalised: GWmodel::gw.dist() returns rp x dp under
  # longlat = TRUE (see .gw_dist_oriented in gw_local_engine.R), which silently
  # corrupted every Great Circle result.
  D <- .gw_dist_oriented(obs_xy, tgt_xy, p, theta, longlat)
  same <- if (!isTRUE(include_self)) outer(obs[[unit]], tgt[[unit]], FUN = "==") else NULL

  Vmat <- as.matrix(as.data.frame(obs)[, value_cols, drop = FALSE])
  storage.mode(Vmat) <- "double"; Vmat[!is.finite(Vmat)] <- 0

  parts <- lapply(kernels, function(km) {
    W <- GWmodel::gw.weight(vdist = D, bw = bw_named[[km]], kernel = km, adaptive = adaptive)
    W <- matrix(W, nrow = nrow(obs_xy), ncol = nrow(tgt_xy))
    if (!is.null(same)) W[same] <- 0
    W[!is.finite(W)] <- 0
    denom <- colSums(W)
    num   <- crossprod(W, Vmat)                            # one BLAS pass, all columns
    bad   <- !is.finite(denom) | denom == 0
    o <- data.table::data.table(uid = tgt[[unit]])
    for (i in seq_along(value_cols)) {
      lag <- num[, i] / denom
      lag[bad] <- NA_real_
      o[[paste0(value_cols[i], "_LM")]] <- lag
    }
    data.table::setnames(o, "uid", unit)
    o[, "kernel" := km]
    o
  })
  out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  data.table::setattr(out, "bandwidth", bw_named)
  out[]
}


#' Geographically weighted spatial lag across kernels (one distance metric)
#'
#' Multi-kernel companion to `estimate_gwlag()`: computes the neighbour-weighted
#' spatial lag under several kernels for a single `distance_metric`, building the
#' (kernel-independent) obs->target distance matrix and the bandwidth only **once**
#' and reusing them across kernels. It is the availability-side analogue of
#' `estimate_gwr_kernels()` / `estimate_gwss_kernels()` - much cheaper than calling
#' `estimate_gwlag()` once per kernel over a spec ensemble.
#'
#' @inheritParams estimate_gwlag
#' @param kernel Character vector of kernels to evaluate (any of `"gaussian"`,
#'   `"exponential"`, `"bisquare"`, `"boxcar"`, `"tricube"`). Default: all five.
#' @param bandwidth Bandwidth strategy when `bw` is `NULL`: `"shared"` (default)
#'   selects one bandwidth (on the first kernel) and reuses it across kernels -
#'   one `GWmodel::bw.gwr()` call instead of `length(kernel)`; `"per_kernel"`
#'   selects a bandwidth per kernel (matches looping `estimate_gwlag()` exactly).
#' @param bw Optional bandwidth override: a single numeric (used for every kernel)
#'   or a **named** numeric keyed by kernel. When supplied no selection is done.
#'
#' @return A `data.table` stacked over `kernel` (a `kernel` column added to the
#'   `estimate_gwlag()` layout: the `unit` id and one `<value_col>_LM` column per
#'   value column). Per-kernel bandwidths are attached as a named numeric in
#'   `attr(., "bandwidth")`.
#' @family Geographically weighted summaries
#' @seealso `estimate_gwlag()`, `estimate_gwr_kernels()`
#' @export
estimate_gwlag_kernels <- function(data, unit, value_cols, geometry = NULL,
                                   coords = c("longitude", "latitude"), poly_id = unit,
                                   predict = NULL, distance_metric = "Euclidean",
                                   kernel = c("gaussian", "exponential", "bisquare",
                                              "boxcar", "tricube"),
                                   adaptive = FALSE, bw = NULL,
                                   bandwidth = c("shared", "per_kernel"),
                                   bw_response = NULL, bw_approach = "CV",
                                   include_self = FALSE) {
  bandwidth <- match.arg(bandwidth)
  dm  <- resolve_distance_metric(distance_metric)
  geo <- .resolve_gw_geometry(data, unit = unit, coords = coords,
                              geometry = geometry, poly_id = poly_id,
                              longlat = dm$longlat, predict = predict)
  .estimate_gwlag_kernels_core(obs = geo$obs, tgt = geo$tgt, unit = unit,
                               value_cols = value_cols, p = dm$p, theta = dm$theta,
                               longlat = dm$longlat, kernels = kernel, adaptive = adaptive,
                               bw = bw, bandwidth = bandwidth, bw_response = bw_response,
                               bw_approach = bw_approach, include_self = include_self)
}
