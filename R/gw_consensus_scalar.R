# ============================================================
# Scalar (continuous) consensus across the distance_metric x kernel domain
# ============================================================
# Sibling of the categorical tools in gw_consensus_class.R, for CONTINUOUS
# per-unit estimates (e.g. a GW coefficient or a warming impact evaluated under
# every kernel x distance_metric setting). Collapses the stack to one robust
# per-unit value with a user-supplied reducer `agg_fun` (default the median),
# and reports the spec-uncertainty spread and a sign-agreement share - the
# scalar analogue of modal agreement. An optional Queen-contiguity smoother
# reconciles each unit's consensus with its neighbourhood.
#
#   consensus       - agg_fun applied across the settings, per unit.
#   sign_agreement  - share of settings whose sign matches the consensus sign.
#   queen_value     - optional: agg_fun over the unit's FIRST-ORDER QUEEN
#                     neighbourhood consensus (spatial de-speckle), widening the
#                     Queen order until at least `min_count` units are gathered
#                     or `max_order` is reached.
#
# Two entry points mirror the GW estimator outputs / the categorical siblings:
#   gw_optimal_scalar_by_point()   - lattice units; Queen from coordinates.
#   gw_optimal_scalar_by_polygon() - polygon units; Queen from shared boundaries.
#
# Shared Queen adjacency builders (.queen_lattice / .queen_polygon) live in
# helpers.R.
# ============================================================


# ------------------------------------------------------------
# Internal: per-unit scalar summary across settings.
# Returns a data.table keyed by .u with the consensus, spread, spec quantiles,
# sign agreement, and the number of settings.
# ------------------------------------------------------------
.scalar_summary <- function(value_dt, unit_col, value_col,
                            agg_fun = stats::median, probs = c(0.05, 0.95)) {
  vd <- data.table::as.data.table(value_dt)
  vd[, .u := as.character(get(unit_col))]
  vd[, .v := as.numeric(get(value_col))]
  vd <- vd[is.finite(.v)]
  vd[, {
    cons <- as.numeric(agg_fun(.v))
    q    <- stats::quantile(.v, probs = probs, names = FALSE, na.rm = TRUE)
    sa   <- if (isTRUE(cons > 0)) mean(.v > 0) else
            if (isTRUE(cons < 0)) mean(.v < 0) else mean(.v == 0)
    list(n_settings     = .N,
         consensus      = cons,
         consensus_sd   = stats::sd(.v),
         consensus_mad  = stats::mad(.v),
         consensus_lo   = q[1L],
         consensus_hi   = q[2L],
         sign_agreement = sa)
  }, by = .u]
}


# ------------------------------------------------------------
# Internal: Queen contiguity smoother for scalars (order expansion).
#   value : per-unit consensus (aligned to `nb`)
#   nb    : list of first-order Queen neighbour indices per unit
# Gathers the unit (optionally) plus its Queen neighbourhood, widening one graph
# hop at a time until at least `min_count` units are gathered or `max_order` is
# reached, then applies `agg_fun` to the gathered values. queen_agreement is the
# neighbourhood share whose sign matches the smoothed value.
# ------------------------------------------------------------
.queen_smooth_scalar <- function(value, nb, agg_fun = stats::median,
                                 include_self = TRUE, min_count = 1L,
                                 max_order = 10L, verbose = FALSE) {
  n   <- length(value)
  val <- rep(NA_real_, n); ord <- rep(NA_integer_, n); agr <- rep(NA_real_, n)
  reduce <- function(cur) {
    v <- value[cur]; v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    qv <- as.numeric(agg_fun(v))
    sa <- if (isTRUE(qv > 0)) mean(v > 0) else
          if (isTRUE(qv < 0)) mean(v < 0) else mean(v == 0)
    list(qv = qv, sa = sa)
  }
  for (i in seq_len(n)) {
    seen     <- if (include_self) i else integer(0)
    frontier <- nb[[i]]
    cur      <- unique(c(seen, frontier))
    order_used <- 1L
    repeat {
      if (!length(cur)) break                          # isolated, no self
      if (length(cur) >= min_count || order_used >= max_order) {
        r <- reduce(cur)
        if (!is.null(r)) { val[i] <- r$qv; agr[i] <- r$sa }
        ord[i] <- order_used; break
      }
      newf <- setdiff(unique(unlist(nb[frontier], use.names = FALSE)), cur)
      if (!length(newf)) {                             # cannot expand further
        r <- reduce(cur)
        if (!is.null(r)) { val[i] <- r$qv; agr[i] <- r$sa }
        ord[i] <- order_used; break
      }
      frontier <- newf; cur <- c(cur, newf); order_used <- order_used + 1L
    }
    if (is.na(val[i])) val[i] <- value[i]              # fallback: keep own value
    if (verbose && i %% 5000L == 0L) message("  queen smooth ", i, " / ", n)
  }
  list(queen_value = val, queen_order = ord, queen_agreement = agr)
}


# ------------------------------------------------------------
# Internal: assemble the scalar output table.
# ------------------------------------------------------------
.assemble_scalar <- function(s, xy, sm, unit_col) {
  out <- data.table::data.table(
    unit_id         = s$.u,
    unit_level      = unit_col,
    longitude       = xy[, 1L],
    latitude        = xy[, 2L],
    n_settings      = s$n_settings,
    consensus       = s$consensus,
    consensus_sd    = round(s$consensus_sd, 6),
    consensus_mad   = round(s$consensus_mad, 6),
    consensus_lo    = s$consensus_lo,
    consensus_hi    = s$consensus_hi,
    sign_agreement  = round(s$sign_agreement, 3)
  )
  if (!is.null(sm)) {
    out[, queen_value     := sm$queen_value]
    out[, queen_order     := sm$queen_order]
    out[, queen_agreement := round(sm$queen_agreement, 3)]
  }
  data.table::setnames(out, "unit_id", unit_col)
  out[]
}


#' Consensus scalar for point (lattice) units
#'
#' Collapses a stacked table of per-unit continuous estimates (one row per unit
#' per `distance_metric` x `kernel` setting, within a history window) to a single
#' per-unit consensus value via a user-supplied reducer (`agg_fun`, default the
#' median), together with the spec-uncertainty spread and a sign-agreement share.
#' Optionally reconciles each unit's consensus with its first-order Queen
#' neighbourhood, widening the Queen order until enough neighbours are gathered.
#' Coordinates are read directly from the table, matching `estimate_gwss_by_point()`.
#'
#' @param value_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`), the continuous estimate (`value_col`), and the two
#'   coordinate columns named in `coords`.
#' @param unit_col Name of the spatial-unit identifier column (e.g. `"grid_id"`).
#' @param value_col Name of the continuous value column. Default `"estimate"`.
#' @param coords Length-2 character vector naming the longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param agg_fun Function reducing a numeric vector (the values across settings
#'   for one unit) to a single scalar. Default `stats::median`. Any function of
#'   the form `function(x) ...` works, e.g. `mean`, `function(x) mean(x, trim = 0.1)`.
#' @param probs Length-2 numeric; the lower/upper quantiles reported across
#'   settings as `consensus_lo`/`consensus_hi`. Default `c(0.05, 0.95)`.
#' @param queen_smooth Logical; if `TRUE`, also compute the Queen-neighbourhood
#'   smoothed consensus. Default `FALSE`.
#' @param include_self Logical; include the unit itself in the Queen
#'   neighbourhood. Default `TRUE`.
#' @param min_count Minimum number of units to gather before reducing the Queen
#'   neighbourhood; the order widens until reached. Default `1`.
#' @param max_order Maximum Queen order to expand to. Default `10`.
#' @param verbose Logical; print progress. Default `FALSE`.
#'
#' @return A `data.table` with one row per unit: the id column, `longitude`,
#'   `latitude`, `n_settings`, `consensus`, `consensus_sd`, `consensus_mad`,
#'   `consensus_lo`, `consensus_hi`, `sign_agreement`, and - when
#'   `queen_smooth = TRUE` - `queen_value`, `queen_order`, `queen_agreement`.
#' @family Consensus scalar
#' @export
gw_optimal_scalar_by_point <- function(value_dt, unit_col,
                                       value_col = "estimate",
                                       coords = c("longitude", "latitude"),
                                       agg_fun = stats::median,
                                       probs = c(0.05, 0.95),
                                       queen_smooth = FALSE,
                                       include_self = TRUE, min_count = 1L,
                                       max_order = 10L, verbose = FALSE) {
  vd <- data.table::as.data.table(value_dt)
  s  <- .scalar_summary(vd, unit_col, value_col, agg_fun, probs)

  co <- unique(vd[, c(unit_col, coords), with = FALSE], by = unit_col)
  data.table::setnames(co, coords, c("longitude", "latitude"))
  co <- co[match(s$.u, as.character(co[[unit_col]]))]
  xy <- as.matrix(co[, .(longitude, latitude)])

  sm <- NULL
  if (isTRUE(queen_smooth)) {
    nb <- .queen_lattice(xy)
    sm <- .queen_smooth_scalar(s$consensus, nb, agg_fun, include_self,
                               min_count, max_order, verbose)
  }
  .assemble_scalar(s, xy, sm, unit_col)
}


#' Consensus scalar for polygon units
#'
#' As `gw_optimal_scalar_by_point()`, but the spatial units are polygons: the
#' Queen neighbours are shared-boundary contiguities of a supplied geometry
#' (matching `estimate_gwss_by_polygon()`, keyed by e.g. `polygon_fips`), and
#' centroids for the output come from the same geometry. Accepts a `terra`
#' `SpatVector` or an `sf` object.
#'
#' @param value_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`) and the continuous estimate (`value_col`).
#' @param unit_col Name of the spatial-unit identifier column (e.g.
#'   `"polygon_fips"`).
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit_col`.
#' @param value_col Name of the continuous value column. Default `"estimate"`.
#' @param poly_id Name of the id field in `polygons`. Default: value of
#'   `unit_col`.
#' @param agg_fun,probs,queen_smooth,include_self,min_count,max_order,verbose See
#'   `gw_optimal_scalar_by_point()`.
#'
#' @return A `data.table` with one row per unit (see
#'   `gw_optimal_scalar_by_point()`).
#' @family Consensus scalar
#' @export
gw_optimal_scalar_by_polygon <- function(value_dt, unit_col, polygons,
                                         value_col = "estimate",
                                         poly_id = unit_col,
                                         agg_fun = stats::median,
                                         probs = c(0.05, 0.95),
                                         queen_smooth = FALSE,
                                         include_self = TRUE, min_count = 1L,
                                         max_order = 10L, verbose = FALSE) {
  s <- .scalar_summary(value_dt, unit_col, value_col, agg_fun, probs)

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
  ci <- match(s$.u, cid)
  xy <- cbind(cxy[ci, 1L], cxy[ci, 2L])

  sm <- NULL
  if (isTRUE(queen_smooth)) {
    nb <- .queen_polygon(polygons, poly_id, s$.u)
    sm <- .queen_smooth_scalar(s$consensus, nb, agg_fun, include_self,
                               min_count, max_order, verbose)
  }
  .assemble_scalar(s, xy, sm, unit_col)
}
