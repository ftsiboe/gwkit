# ============================================================
# Consensus across the distance_metric x kernel domain
# ============================================================
# Reduces results estimated across many settings (every kernel x distance_metric
# of a GW estimator, within one history window) to a single per-unit summary.
# Two flavours share this file:
#
#   gw_consensus_scalar() - CONTINUOUS estimates (e.g. a GW coefficient / impact):
#       a robust per-unit consensus via a user-supplied reducer `agg_fun`
#       (default the median), the spec-uncertainty spread, a sign-agreement share,
#       and an optional Queen-contiguity smoother.
#   gw_consensus_class()  - CATEGORICAL classifications: the modal class across
#       settings plus a Queen-contiguity vote (tie-driven order expansion) for
#       optional spatial de-speckling.
#
# Point vs polygon mode is class-detected in both (an sf/SpatVector `geometry` ->
# polygon units; else `coords` columns of the table). The shared point/Queen
# resolution is .consensus_geometry(); the Queen adjacency builders
# (.queen_lattice / .queen_polygon) and the representative-point helper
# (.gwkit_centroids) live in helpers.R / estimate_gwlag.R.
# ============================================================

utils::globalVariables(c(".u", ".c", ".v", ".N", "longitude", "latitude",
                         "queen_value", "queen_order", "queen_agreement"))


# ------------------------------------------------------------
# Internal: shared point + Queen resolution for the consensus estimators.
# Class-detected: a polygon `geometry` gives point-on-surface coordinates (from
# .gwkit_centroids) and shared-boundary Queen adjacency; otherwise the `coords`
# columns of `dt` give coordinates and lattice Queen adjacency. `queen_smooth`
# controls whether the neighbour list `nb` is built. Returns list(xy, nb),
# aligned to `unit_ids`.
# ------------------------------------------------------------
.consensus_geometry <- function(unit_ids, dt, unit_col, geometry, coords,
                                poly_id, queen_smooth) {
  is_poly <- inherits(geometry, c("sf", "SpatVector"))
  if (is_poly) {
    cds <- .gwkit_centroids(geometry, poly_id = poly_id, longlat = TRUE)
    ci  <- match(unit_ids, as.character(cds[[".gw_uid"]]))
    xy  <- cbind(cds[[".gw_x"]][ci], cds[[".gw_y"]][ci])
  } else if (!is.null(geometry)) {
    stop("`geometry` must be a terra SpatVector or an sf object.")
  } else {
    co <- unique(dt[, c(unit_col, coords), with = FALSE], by = unit_col)
    data.table::setnames(co, coords, c("longitude", "latitude"))
    co <- co[match(unit_ids, as.character(co[[unit_col]]))]
    xy <- as.matrix(co[, .(longitude, latitude)])
  }
  nb <- if (isTRUE(queen_smooth)) {
    if (is_poly) .queen_polygon(geometry, poly_id, unit_ids) else .queen_lattice(xy)
  } else NULL
  list(xy = xy, nb = nb)
}


# ============================================================
# Scalar (continuous) consensus
# ============================================================

# ------------------------------------------------------------
# Internal: per-unit scalar summary across settings. Returns a data.table keyed
# by .u with the consensus, spread, spec quantiles, sign agreement, n settings.
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


#' Consensus scalar across the kernel x distance-metric domain
#'
#' Collapses a stacked table of per-unit continuous estimates (one row per unit
#' per `distance_metric` x `kernel` setting, within a history window) to a single
#' per-unit consensus value via a user-supplied reducer (`agg_fun`, default the
#' median), together with the spec-uncertainty spread and a sign-agreement share.
#' Optionally reconciles each unit's consensus with its first-order Queen
#' neighbourhood, widening the Queen order until enough neighbours are gathered.
#' Point vs polygon mode is class-detected: pass an `sf`/`SpatVector` `geometry`
#' for polygon units (Queen = shared boundaries, output coords = point-on-surface);
#' otherwise the `coords` columns of `value_dt` are used (Queen = lattice).
#'
#' @param value_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`), the continuous estimate (`value_col`), and - in
#'   point mode - the coordinate columns named in `coords`.
#' @param unit_col Name of the spatial-unit identifier column (e.g. `"grid_id"`).
#' @param geometry Optional `sf`/`SpatVector` polygon layer (id field `poly_id`
#'   matching `unit_col`) for polygon mode. If `NULL`, point mode uses `coords`.
#' @param value_col Name of the continuous value column. Default `"estimate"`.
#' @param coords Length-2 character vector naming the longitude/latitude columns
#'   (point mode). Default `c("longitude", "latitude")`.
#' @param poly_id Name of the id field in `geometry`. Default: `unit_col`.
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
#' @family Consensus
#' @export
gw_consensus_scalar <- function(value_dt, unit_col, geometry = NULL,
                              value_col = "estimate",
                              coords = c("longitude", "latitude"),
                              poly_id = unit_col,
                              agg_fun = stats::median, probs = c(0.05, 0.95),
                              queen_smooth = FALSE, include_self = TRUE,
                              min_count = 1L, max_order = 10L, verbose = FALSE) {
  vd  <- data.table::as.data.table(value_dt)
  s   <- .scalar_summary(vd, unit_col, value_col, agg_fun, probs)
  geo <- .consensus_geometry(s$.u, vd, unit_col, geometry, coords, poly_id,
                             queen_smooth = queen_smooth)
  sm  <- if (isTRUE(queen_smooth))
           .queen_smooth_scalar(s$consensus, geo$nb, agg_fun, include_self,
                                min_count, max_order, verbose) else NULL
  .assemble_scalar(s, geo$xy, sm, unit_col)
}


# ============================================================
# Categorical (class) consensus
# ============================================================

# ------------------------------------------------------------
# Internal: per-unit class counts across settings -> modal class.
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
# Internal: assemble the class output table.
# ------------------------------------------------------------
.assemble_consensus <- function(cc, xy, vote, unit_col) {
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


#' Consensus class across the kernel x distance-metric domain
#'
#' Collapses a stacked table of per-unit classifications (one row per unit per
#' `distance_metric` x `kernel` setting, within a history window) to a single
#' consensus class per unit: the mode across settings (`modal_class`) plus a
#' first-order Queen contiguity vote for optional spatial de-speckling
#' (`queen_class`), with ties resolved by widening the Queen order. Point vs
#' polygon mode is class-detected: pass an `sf`/`SpatVector` `geometry` for
#' polygon units (Queen = shared boundaries); otherwise `coords` columns are read
#' from the table (Queen = lattice).
#'
#' @param class_dt A `data.table`/data frame stacked over settings, containing
#'   the unit id (`unit_col`), the categorical class (`class_col`), and - in
#'   point mode - the coordinate columns named in `coords`.
#' @param unit_col Name of the spatial-unit identifier column (e.g. `"grid_id"`).
#' @param geometry Optional `sf`/`SpatVector` polygon layer (id field `poly_id`
#'   matching `unit_col`) for polygon mode. If `NULL`, point mode uses `coords`.
#' @param class_col Name of the categorical class column. Default
#'   `"trend_class_05"`.
#' @param coords Length-2 character vector naming the longitude/latitude columns
#'   (point mode). Default `c("longitude", "latitude")`.
#' @param poly_id Name of the id field in `geometry`. Default: `unit_col`.
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
#' @family Consensus
#' @export
gw_consensus_class <- function(class_dt, unit_col, geometry = NULL,
                             class_col = "trend_class_05",
                             coords = c("longitude", "latitude"),
                             poly_id = unit_col,
                             include_self = TRUE, max_order = 10L,
                             class_levels = NULL, verbose = FALSE) {
  cd   <- data.table::as.data.table(class_dt)
  cc   <- .class_counts(cd, unit_col, class_col, class_levels)
  geo  <- .consensus_geometry(cc$units, cd, unit_col, geometry, coords, poly_id,
                              queen_smooth = TRUE)
  vote <- .queen_vote(cc$modal_code, geo$nb, length(cc$class_levels), cc$class_levels,
                      include_self = include_self, max_order = max_order,
                      verbose = verbose)
  .assemble_consensus(cc, geo$xy, vote, unit_col)
}
