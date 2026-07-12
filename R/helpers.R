# ============================================================
# Shared spatial + adjacency helpers
# ============================================================
# Geometry primitives shared across the gwkit estimators and consensus tools:
#
#   .gwkit_centroids()     - polygon -> representative point (point-on-surface).
#   .resolve_gw_geometry() - class-detected point vs polygon geometry for the
#                            regression estimators; returns obs/tgt coordinates.
#   .queen_lattice()       - first-order Queen adjacency on a regular lattice.
#   .queen_polygon()       - first-order Queen adjacency for polygons (shared
#                            boundary contiguity).
# ============================================================

utils::globalVariables(c(".u", ".c", ".v", ".N", "longitude", "latitude"))

# Internal coordinate/id column names used to shuttle geometry between the
# resolver and the estimator cores. Dot-prefixed so they never collide with a
# user's model columns (a variable literally named `X`/`Y`/`Xo`... is legal).
#   .gw_uid : unit id from the polygon layer
#   .gw_x / .gw_y   : representative-point coordinates (centroids output)
#   .gw_xo / .gw_yo : source (obs) coordinates ; .gw_xt / .gw_yt : target coords


# ------------------------------------------------------------
# Internal: representative point coordinates (X, Y) for polygon units, projected
# to EPSG:4326 when a longlat metric is requested. Returns data.table(uid, X, Y).
# The representative point is the POINT-ON-SURFACE (guaranteed inside the polygon)
# - terra::centroids(inside = TRUE) / sf::st_point_on_surface() - so every gwkit
# polygon estimator reduces a polygon to the same point.
# ------------------------------------------------------------
.gwkit_centroids <- function(polygons, poly_id, longlat) {
  if (inherits(polygons, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required for SpatVector polygons.")
    pg  <- if (isTRUE(longlat) && !terra::is.lonlat(polygons))
             terra::project(polygons, "EPSG:4326") else polygons
    pt  <- terra::centroids(pg, inside = TRUE); xy <- terra::crds(pt)
    data.table::data.table(.gw_uid = as.character(terra::values(pt)[[poly_id]]),
                           .gw_x = xy[, 1L], .gw_y = xy[, 2L])
  } else if (inherits(polygons, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for sf polygons.")
    pg  <- if (isTRUE(longlat)) sf::st_transform(polygons, 4326) else polygons
    pt  <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(pg)))
    xy  <- sf::st_coordinates(pt)
    data.table::data.table(.gw_uid = as.character(sf::st_drop_geometry(pg)[[poly_id]]),
                           .gw_x = xy[, 1L], .gw_y = xy[, 2L])
  } else {
    stop("`polygons` must be a terra SpatVector or an sf object.")
  }
}


# ------------------------------------------------------------
# Internal: class-detected geometry resolver for the regression estimators. A
# polygon layer (from `geometry`, or from `data` itself when it is sf/SpatVector)
# reduces each unit to its point-on-surface via .gwkit_centroids(); otherwise the
# `coords` columns of `data` are used directly. Returns obs (unit, .gw_xo,
# .gw_yo, + model cols) and tgt (unit, .gw_xt, .gw_yt). `predict` restricts: a
# character/id vector (or a data frame with `unit`) in polygon mode, a data frame
# of point targets in point mode; NULL uses the unique units in `data`.
# ------------------------------------------------------------
.resolve_gw_geometry <- function(data, unit, coords = c("longitude", "latitude"),
                                 geometry = NULL, poly_id = unit, longlat = FALSE,
                                 predict = NULL) {
  if (!is.null(geometry) && !inherits(geometry, c("sf", "SpatVector")))
    stop("`geometry` must be a terra SpatVector or an sf object.")
  poly <- NULL
  if (inherits(geometry, c("sf", "SpatVector")))    poly <- geometry
  else if (inherits(data, c("sf", "SpatVector")))   poly <- data

  if (!is.null(poly)) {
    cds <- .gwkit_centroids(poly, poly_id = poly_id, longlat = longlat)  # .gw_uid/.gw_x/.gw_y
    d <- if (inherits(data, "sf")) {
           data.table::as.data.table(sf::st_drop_geometry(data))
         } else if (inherits(data, "SpatVector")) {
           data.table::as.data.table(terra::values(data))
         } else data.table::as.data.table(data)
    d[[unit]] <- as.character(d[[unit]])
    n_units <- data.table::uniqueN(d[[unit]])
    obs <- merge(d, cds, by.x = unit, by.y = ".gw_uid", all.x = FALSE)
    n_dropped <- n_units - data.table::uniqueN(obs[[unit]])
    if (n_dropped > 0L)
      message(sprintf(
        "gwkit: %d of %d unit(s) had no matching polygon in `geometry` (poly_id = '%s') and were dropped.",
        n_dropped, n_units, poly_id))
    data.table::setnames(obs, c(".gw_x", ".gw_y"), c(".gw_xo", ".gw_yo"))
    tgt <- data.table::copy(cds)
    data.table::setnames(tgt, c(".gw_uid", ".gw_x", ".gw_y"), c(unit, ".gw_xt", ".gw_yt"))
    if (!is.null(predict)) {
      pid <- if (is.data.frame(predict)) as.character(predict[[unit]]) else as.character(predict)
      tgt <- tgt[get(unit) %in% pid]
    }
  } else {
    d <- data.table::as.data.table(data); d[[unit]] <- as.character(d[[unit]])
    obs <- data.table::copy(d)
    obs[, c(".gw_xo", ".gw_yo") := list(as.numeric(d[[coords[1L]]]),
                                        as.numeric(d[[coords[2L]]]))]
    if (is.null(predict)) {
      tgt <- unique(data.table::copy(d), by = unit)
      tgt[, c(".gw_xt", ".gw_yt") := list(as.numeric(tgt[[coords[1L]]]),
                                          as.numeric(tgt[[coords[2L]]]))]
    } else {
      pd <- data.table::as.data.table(predict); pd[[unit]] <- as.character(pd[[unit]])
      tgt <- pd
      tgt[, c(".gw_xt", ".gw_yt") := list(as.numeric(pd[[coords[1L]]]),
                                          as.numeric(pd[[coords[2L]]]))]
    }
  }
  list(obs = obs, tgt = tgt, mode = if (!is.null(poly)) "polygon" else "point")
}


# ------------------------------------------------------------
# Internal: first-order Queen adjacency on a regular lattice, from centroids.
# Two cells are Queen neighbours if they are within one step in both lon and
# lat (the 8 surrounding cells). Step sizes are inferred from the coordinates.
# ------------------------------------------------------------
.queen_lattice <- function(xy) {
  lon <- xy[, 1L]; lat <- xy[, 2L]
  ulon <- sort(unique(round(lon, 6))); ulat <- sort(unique(round(lat, 6)))
  step_lon <- if (length(ulon) > 1L) min(diff(ulon)) else 1
  step_lat <- if (length(ulat) > 1L) min(diff(ulat)) else 1
  ix <- as.integer(round((lon - min(lon, na.rm = TRUE)) / step_lon))
  iy <- as.integer(round((lat - min(lat, na.rm = TRUE)) / step_lat))
  key <- paste(ix, iy, sep = "_")
  lut <- stats::setNames(seq_along(key), key)
  offs <- as.matrix(expand.grid(dx = -1:1, dy = -1:1))
  offs <- offs[!(offs[, 1L] == 0L & offs[, 2L] == 0L), , drop = FALSE]
  nb <- vector("list", length(key))
  for (m in seq_along(key)) {
    if (is.na(ix[m]) || is.na(iy[m])) { nb[[m]] <- integer(0L); next }
    nk  <- paste(ix[m] + offs[, 1L], iy[m] + offs[, 2L], sep = "_")
    hit <- lut[nk]
    nb[[m]] <- unname(hit[!is.na(hit)])
  }
  nb
}


# ------------------------------------------------------------
# Internal: first-order Queen adjacency for polygons (shared boundary),
# returned aligned to `units`. Accepts terra SpatVector or sf.
# ------------------------------------------------------------
.queen_polygon <- function(polygons, poly_id, units) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for polygon Queen contiguity.")
  pg <- if (inherits(polygons, "SpatVector")) sf::st_as_sf(polygons) else polygons
  ids <- as.character(sf::st_drop_geometry(pg)[[poly_id]])
  touch <- sf::st_touches(pg)                          # first-order Queen (sparse)
  pos <- match(units, ids)                             # polygon row for each unit
  nb <- vector("list", length(units))
  for (m in seq_along(units)) {
    pj <- pos[m]
    if (is.na(pj)) { nb[[m]] <- integer(0L); next }
    nid <- ids[touch[[pj]]]
    idx <- match(nid, units)
    nb[[m]] <- idx[!is.na(idx)]
  }
  nb
}
