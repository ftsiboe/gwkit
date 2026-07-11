# ============================================================
# Optimal-class selection across the distance_metric x kernel domain
# ============================================================
# Given a stacked table of per-unit categorical classifications produced across
# many settings (every kernel x distance_metric of a GW estimator, within one
# history window), pick a single "optimal" class per spatial unit:
#
#   modal_class - the mode across the distance_metric x kernel settings (per
#                 unit). This IS the optimal class: the consensus pick over the
#                 setting domain.
#   queen_class - a contiguity vote for optional spatial de-speckling: the
#                 majority modal class among a unit's FIRST-ORDER QUEEN
#                 neighbours; when that majority ties, the neighbourhood is
#                 widened one order at a time (2nd, 3rd, ...) until the tie
#                 breaks or `max_order` is reached.
#
# Everything is agnostic to the spatial-unit identifier: pass the id column via
# `unit_col` and the class column via `class_col`. Two entry points match the GW
# estimator outputs:
#
#   gw_optimal_class_by_point()   - units on a regular lattice; Queen adjacency
#                                   is derived from the point coordinates (e.g.
#                                   the output of estimate_gwss_by_point()).
#   gw_optimal_class_by_polygon() - polygon units; Queen adjacency is the shared-
#                                   boundary contiguity of the geometry (e.g.
#                                   estimate_gwss_by_polygon() keyed by
#                                   polygon_fips).
#
# Shared Queen adjacency builders (.queen_lattice / .queen_polygon) live in
# helpers.R.
# ============================================================


# ------------------------------------------------------------
# Internal: per-unit class counts across settings -> modal class.
# Returns units (ordered), the count matrix N, modal code/label/agreement,
# and the class levels used.
# ------------------------------------------------------------
.class_counts <- function(class_dt, unit_col, class_col, class_levels = NULL) {
  cd <- data.table::as.data.table(class_dt)
  cd <- cd[!is.na(get(class_col))]
  cd[, .u := as.character(get(unit_col))]
  if (is.null(class_levels))
    class_levels <- sort(unique(as.character(cd[[class_col]])))
  cd[, .c := factor(as.character(get(class_col)), levels = class_levels)]

  counts <- cd[, .N, by = c(".u", ".c")]
  wide   <- data.table::dcast(counts, .u ~ .c, value.var = "N",
                              fill = 0, drop = c(TRUE, FALSE))
  units  <- wide$.u
  N      <- as.matrix(wide[, ..class_levels]); storage.mode(N) <- "double"
  row_max <- N[cbind(seq_len(nrow(N)), max.col(N, ties.method = "first"))]

  list(units = units, N = N, class_levels = class_levels,
       modal_code  = max.col(N, ties.method = "first"),
       modal_class = class_levels[max.col(N, ties.method = "first")],
       modal_agree = row_max / pmax(rowSums(N), 1))
}


# ------------------------------------------------------------
# Internal: Queen contiguity vote with tie-driven order expansion.
#   modal_code : integer class code per unit (aligned to `nb`)
#   nb         : list of first-order Queen neighbour indices per unit
# Starts at first-order Queen; while the winning class ties, widens the
# neighbourhood one graph hop at a time until the tie breaks or max_order.
# ------------------------------------------------------------
.queen_vote <- function(modal_code, nb, nlev, class_levels,
                        include_self = TRUE, max_order = 10L, verbose = FALSE) {
  n <- length(modal_code)
  cls <- rep(NA_character_, n); ord <- rep(NA_integer_, n); agr <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    seen     <- if (include_self) i else integer(0)
    frontier <- nb[[i]]
    cur      <- unique(c(seen, frontier))
    order_used <- 1L
    repeat {
      if (!length(cur)) { break }                      # isolated, no self
      tab     <- tabulate(modal_code[cur], nbins = nlev)
      mx      <- max(tab); winners <- which(tab == mx)
      resolved <- length(winners) == 1L
      if (resolved || order_used >= max_order) {
        cls[i] <- class_levels[winners[1L]]            # first level breaks a final tie
        ord[i] <- order_used
        agr[i] <- mx / length(cur)
        break
      }
      # widen one Queen order (graph hop) and retry
      newf <- setdiff(unique(unlist(nb[frontier], use.names = FALSE)), cur)
      if (!length(newf)) {                             # cannot expand further
        cls[i] <- class_levels[winners[1L]]; ord[i] <- order_used
        agr[i] <- mx / length(cur); break
      }
      frontier <- newf; cur <- c(cur, newf); order_used <- order_used + 1L
    }
    if (verbose && i %% 5000L == 0L) message("  queen vote ", i, " / ", n)
  }
  list(queen_class = cls, queen_order = ord, queen_agreement = agr)
}


# ------------------------------------------------------------
# Internal: assemble the output table.
# ------------------------------------------------------------
.assemble_optimal <- function(cc, xy, vote, unit_col) {
  out <- data.table::data.table(
    unit_id         = cc$units,
    unit_level      = unit_col,
    longitude       = xy[, 1L],
    latitude        = xy[, 2L],
    n_settings      = rowSums(cc$N),
    modal_class     = cc$modal_class,
    modal_agreement = round(cc$modal_agree, 3),
    queen_class     = vote$queen_class,
    queen_order     = vote$queen_order,
    queen_agreement = round(vote$queen_agreement, 3)
  )
  data.table::setnames(out, "unit_id", unit_col)
  out[]
}


#' Optimal class for point (lattice) units
#'
#' Collapses a stacked table of per-unit classifications (one row per unit per
#' `distance_metric` x `kernel` setting, within a history window) to a single
#' optimal class per unit: the mode across settings (`modal_class`) plus a
#' first-order Queen contiguity vote for optional spatial de-speckling
#' (`queen_class`), with ties resolved by widening the Queen order. Coordinates
#' are read directly from the table, matching `estimate_gwss_by_point()`.
#'
#' @param class_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`), the categorical class (`class_col`), and the two
#'   coordinate columns named in `coords`.
#' @param unit_col Name of the spatial-unit identifier column (e.g. `"grid_id"`).
#' @param class_col Name of the categorical class column. Default
#'   `"trend_class_05"`.
#' @param coords Length-2 character vector naming the longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param include_self Logical; count the unit itself in the Queen vote. Default
#'   `TRUE`.
#' @param max_order Maximum Queen order to expand to when breaking ties. Default
#'   `10`.
#' @param class_levels Optional character vector of all class labels (fixes the
#'   count columns / vote ordering). Default: sorted unique observed classes.
#' @param verbose Logical; print progress. Default `FALSE`.
#'
#' @return A `data.table` with one row per unit: the id column, `longitude`,
#'   `latitude`, `n_settings`, `modal_class`, `modal_agreement`, `queen_class`,
#'   `queen_order` (order at which the vote resolved), and `queen_agreement`.
#' @family Optimal class selection
#' @export
gw_optimal_class_by_point <- function(class_dt, unit_col,
                                      class_col = "trend_class_05",
                                      coords = c("longitude", "latitude"),
                                      include_self = TRUE, max_order = 10L,
                                      class_levels = NULL, verbose = FALSE) {
  cd <- data.table::as.data.table(class_dt)
  cc <- .class_counts(cd, unit_col, class_col, class_levels)

  co <- unique(cd[, c(unit_col, coords), with = FALSE], by = unit_col)
  data.table::setnames(co, coords, c("longitude", "latitude"))
  co <- co[match(cc$units, as.character(co[[unit_col]]))]
  xy <- as.matrix(co[, .(longitude, latitude)])

  nb   <- .queen_lattice(xy)
  vote <- .queen_vote(cc$modal_code, nb, length(cc$class_levels), cc$class_levels,
                      include_self = include_self, max_order = max_order,
                      verbose = verbose)
  .assemble_optimal(cc, xy, vote, unit_col)
}


#' Optimal class for polygon units
#'
#' As `gw_optimal_class_by_point()`, but the spatial units are polygons: the
#' Queen neighbours are shared-boundary contiguities of a supplied geometry
#' (matching `estimate_gwss_by_polygon()`, keyed by e.g. `polygon_fips`), and
#' centroids for the output come from the same geometry. Accepts a `terra`
#' `SpatVector` or an `sf` object.
#'
#' @param class_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`) and the categorical class (`class_col`).
#' @param unit_col Name of the spatial-unit identifier column (e.g.
#'   `"polygon_fips"`).
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit_col`.
#' @param class_col Name of the categorical class column. Default
#'   `"trend_class_05"`.
#' @param poly_id Name of the id field in `polygons`. Default: value of
#'   `unit_col`.
#' @param include_self,max_order,class_levels,verbose See
#'   `gw_optimal_class_by_point()`.
#'
#' @return A `data.table` with one row per unit (see
#'   `gw_optimal_class_by_point()`).
#' @family Optimal class selection
#' @export
gw_optimal_class_by_polygon <- function(class_dt, unit_col, polygons,
                                        class_col = "trend_class_05",
                                        poly_id = unit_col,
                                        include_self = TRUE, max_order = 10L,
                                        class_levels = NULL, verbose = FALSE) {
  cc <- .class_counts(class_dt, unit_col, class_col, class_levels)

  # centroids (for output coordinates), aligned to units
  if (inherits(polygons, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required for SpatVector polygons.")
    pg  <- if (!terra::is.lonlat(polygons)) terra::project(polygons, "EPSG:4326") else polygons
    cen <- terra::centroids(pg); cxy <- terra::crds(cen)
    cid <- as.character(terra::values(cen)[[poly_id]])
  } else if (inherits(polygons, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for sf polygons.")
    pg  <- sf::st_transform(polygons, 4326)
    cen <- suppressWarnings(sf::st_centroid(sf::st_geometry(pg)))
    cxy <- sf::st_coordinates(cen)
    cid <- as.character(sf::st_drop_geometry(pg)[[poly_id]])
  } else {
    stop("`polygons` must be a terra SpatVector or an sf object.")
  }
  ci <- match(cc$units, cid)
  xy <- cbind(cxy[ci, 1L], cxy[ci, 2L])

  nb   <- .queen_polygon(polygons, poly_id, cc$units)
  vote <- .queen_vote(cc$modal_code, nb, length(cc$class_levels), cc$class_levels,
                      include_self = include_self, max_order = max_order,
                      verbose = verbose)
  .assemble_optimal(cc, xy, vote, unit_col)
}
