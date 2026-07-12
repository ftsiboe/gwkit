# ============================================================
# Geographically weighted spatial lag (neighbour-weighted mean)
# ============================================================
# estimate_gwlag_by_point()   - units carry point coordinates.
# estimate_gwlag_by_polygon() - units are polygons; centroids from geometry.
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

utils::globalVariables(c("Xo", "Yo", "Xt", "Yt", "uid", "X", "Y"))

# ------------------------------------------------------------
# Internal: centroid coordinates (x, y) for polygon units, projected to
# EPSG:4326 when a longlat metric is requested. Returns data.table(uid, X, Y).
# ------------------------------------------------------------
.gwkit_centroids <- function(polygons, poly_id, longlat) {
  if (inherits(polygons, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required for SpatVector polygons.")
    pg  <- if (isTRUE(longlat) && !terra::is.lonlat(polygons))
             terra::project(polygons, "EPSG:4326") else polygons
    cen <- terra::centroids(pg); xy <- terra::crds(cen)
    data.table::data.table(uid = as.character(terra::values(cen)[[poly_id]]),
                           X = xy[, 1L], Y = xy[, 2L])
  } else if (inherits(polygons, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for sf polygons.")
    pg  <- if (isTRUE(longlat)) sf::st_transform(polygons, 4326) else polygons
    cen <- suppressWarnings(sf::st_centroid(sf::st_geometry(pg)))
    xy  <- sf::st_coordinates(cen)
    data.table::data.table(uid = as.character(sf::st_drop_geometry(pg)[[poly_id]]),
                           X = xy[, 1L], Y = xy[, 2L])
  } else {
    stop("`polygons` must be a terra SpatVector or an sf object.")
  }
}


# ------------------------------------------------------------
# Internal engine. `obs` carries the source coordinates (Xo, Yo) and value_cols;
# `tgt` carries the target coordinates (Xt, Yt) and the target ids.
# ------------------------------------------------------------
.estimate_gwlag_core <- function(obs, tgt, unit, value_cols,
                                 p, theta, longlat, kernel, adaptive, bw,
                                 bw_response, bw_approach, include_self) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwlag_*().")

  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)

  obs_xy <- as.matrix(obs[, c("Xo", "Yo")])
  tgt_xy <- as.matrix(tgt[, c("Xt", "Yt")])

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


#' Geographically weighted spatial lag (neighbour-weighted mean) at points
#'
#' For each target unit, computes the geographically weighted mean of one or
#' more value columns over the source units, using `GWmodel::gw.weight()` weights
#' on a `GWmodel::gw.dist()` distance matrix. With `include_self = FALSE` the
#' focal unit is excluded from its own neighbourhood, yielding a strict
#' neighbour mean (a spatial lag).
#'
#' @param data A `data.table`/data frame of the SOURCE units, one row per unit,
#'   containing `unit`, the coordinate columns in `coords`, and all `value_cols`.
#' @param unit Character; the unit-id column (present in both source and target).
#' @param value_cols Character vector of numeric columns to lag. Each returns a
#'   column named `<value_col>_LM`.
#' @param coords Length-2 character vector naming longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param predict_data Optional TARGET units (one row per unit) with `unit` and
#'   `coords`. Defaults to `data` (targets = sources).
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `FALSE`.
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
estimate_gwlag_by_point <- function(data, unit, value_cols,
                                    coords = c("longitude", "latitude"),
                                    predict_data = NULL,
                                    distance_metric = "Euclidean",
                                    kernel = "bisquare", adaptive = FALSE,
                                    bw = NULL, bw_response = NULL,
                                    bw_approach = "CV", include_self = FALSE) {
  dm <- resolve_distance_metric(distance_metric)
  d  <- data.table::as.data.table(data)
  d[[unit]] <- as.character(d[[unit]])
  obs <- data.table::data.table(d)
  obs[, `:=`(Xo = as.numeric(d[[coords[1L]]]), Yo = as.numeric(d[[coords[2L]]]))]

  if (is.null(predict_data)) {
    tgt <- data.table::data.table(d)
    tgt[, `:=`(Xt = as.numeric(d[[coords[1L]]]), Yt = as.numeric(d[[coords[2L]]]))]
  } else {
    pd <- data.table::as.data.table(predict_data)
    pd[[unit]] <- as.character(pd[[unit]])
    tgt <- data.table::data.table(pd)
    tgt[, `:=`(Xt = as.numeric(pd[[coords[1L]]]), Yt = as.numeric(pd[[coords[2L]]]))]
  }

  .estimate_gwlag_core(obs = obs, tgt = tgt, unit = unit, value_cols = value_cols,
                       p = dm$p, theta = dm$theta, longlat = dm$longlat,
                       kernel = kernel, adaptive = adaptive, bw = bw,
                       bw_response = bw_response, bw_approach = bw_approach,
                       include_self = include_self)
}


#' Geographically weighted spatial lag on polygons
#'
#' As `estimate_gwlag_by_point()`, but the spatial units are polygons; centroid
#' coordinates are derived from `polygons` and joined to `data` (and, if given,
#' `predict_data`) by the unit id.
#'
#' @param data Source units (one row per unit) with `unit` and all `value_cols`.
#' @param unit Character; the unit-id column.
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit`. Provides target ids and centroids; the value
#'   columns are joined from `data` (and `predict_ids` restricts the targets).
#' @param value_cols Character vector of numeric columns to lag.
#' @param poly_id Name of the id field in `polygons`. Default: value of `unit`.
#' @param predict_ids Optional character vector of target unit ids (subset of the
#'   polygon ids). Default: all polygon units.
#' @param distance_metric,kernel,adaptive,bw,bw_response,bw_approach,include_self
#'   See `estimate_gwlag_by_point()`.
#'
#' @return A `data.table` with one row per target unit (see
#'   `estimate_gwlag_by_point()`).
#' @family Geographically weighted summaries
#' @export
estimate_gwlag_by_polygon <- function(data, unit, polygons, value_cols,
                                      poly_id = unit, predict_ids = NULL,
                                      distance_metric = "Euclidean",
                                      kernel = "bisquare", adaptive = FALSE,
                                      bw = NULL, bw_response = NULL,
                                      bw_approach = "CV", include_self = FALSE) {
  dm  <- resolve_distance_metric(distance_metric)
  cds <- .gwkit_centroids(polygons, poly_id = poly_id, longlat = dm$longlat)

  d <- data.table::as.data.table(data)
  d[[unit]] <- as.character(d[[unit]])
  # source: units that carry data, joined to their centroid
  obs <- merge(d, cds, by.x = unit, by.y = "uid", all.x = FALSE)
  data.table::setnames(obs, c("X", "Y"), c("Xo", "Yo"))

  # target: all polygon units (or a supplied subset)
  tgt <- data.table::copy(cds)
  data.table::setnames(tgt, c("uid", "X", "Y"), c(unit, "Xt", "Yt"))
  if (!is.null(predict_ids)) tgt <- tgt[get(unit) %in% as.character(predict_ids)]

  .estimate_gwlag_core(obs = obs, tgt = tgt, unit = unit, value_cols = value_cols,
                       p = dm$p, theta = dm$theta, longlat = dm$longlat,
                       kernel = kernel, adaptive = adaptive, bw = bw,
                       bw_response = bw_response, bw_approach = bw_approach,
                       include_self = include_self)
}
