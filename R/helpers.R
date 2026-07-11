# ============================================================
# Consensus across the distance_metric x kernel domain - shared helpers
# ============================================================
# The consensus layer reduces results estimated across many settings (every
# kernel x distance_metric of a GW estimator, within one history window) to a
# single per-unit summary. Two flavours are provided in sibling files:
#
#   gw_consensus_class.R  - categorical outcomes (modal class + Queen vote).
#   gw_consensus_scalar.R - continuous outcomes (a reducer such as the median,
#                           spec-uncertainty spread, sign agreement, and an
#                           optional Queen smoother).
#
# Both flavours share the first-order Queen adjacency builders defined here:
#
#   .queen_lattice()  - Queen adjacency on a regular lattice, from centroids.
#   .queen_polygon()  - Queen adjacency for polygons (shared-boundary contiguity).
# ============================================================

utils::globalVariables(c(".u", ".c", ".v", ".N", "longitude", "latitude"))


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
