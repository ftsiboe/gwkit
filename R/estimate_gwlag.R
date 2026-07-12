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
  D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy,
                        focus = 0, p = p, theta = theta, longlat = longlat)
  W <- GWmodel::gw.weight(vdist = D, bw = bw, kernel = kernel, adaptive = adaptive)
  W <- matrix(W, nrow = nrow(obs_xy), ncol = nrow(tgt_xy))   # obs x target

  # strict neighbour mean: zero a source's weight onto its own target
  if (!isTRUE(include_self)) {
    same <- outer(obs[[unit]], tgt[[unit]], FUN = "==")
    W[same] <- 0
  }

  denom <- colSums(W, na.rm = TRUE)
  out <- data.table::data.table(uid = tgt[[unit]])
  for (vc in value_cols) {
    z <- as.numeric(obs[[vc]]); z[!is.finite(z)] <- 0
    num <- colSums(W * z, na.rm = TRUE)
    lag <- num / denom
    lag[!is.finite(denom) | denom == 0] <- NA_real_
    out[[paste0(vc, "_LM")]] <- lag
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
