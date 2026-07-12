# ============================================================
# Geographically weighted local regression coefficients (multi-covariate)
# ============================================================
# estimate_gwr_coefficients_by_point()   - units carry point coordinates.
# estimate_gwr_coefficients_by_polygon() - units are polygons; centroids derived.
#
# For every target unit i, fits a locally weighted linear model
#   lm(formula, weights = w_i)
# where w_i are GWmodel::gw.weight() weights over the source units, and returns
# the full coefficient table (one row per model term). Unlike
# estimate_gwr_by_point() (a single-covariate mean/variance slope with an
# analytic SE), this accepts an ARBITRARY formula (e.g. y ~ x1 + x2) and returns
# exactly what stats::lm(weights=) reports - so it reproduces a hand-rolled
# per-unit weighted lm() bit for bit, including standard errors and residual df.
# ============================================================

utils::globalVariables(c("Xo", "Yo", "Xt", "Yt", "uid", "X", "Y", ".w",
                         "term", "unit_level", "model_estimator"))

# ------------------------------------------------------------
# Internal engine. `obs` carries source coordinates (Xo, Yo) and the model
# variables; `tgt` carries target coordinates (Xt, Yt) and target ids.
# ------------------------------------------------------------
.estimate_gwrcoef_core <- function(obs, tgt, unit, formula,
                                   p, theta, longlat, kernel, adaptive, bw,
                                   bw_response, bw_approach, terms) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwr_coefficients_*().")

  formula <- stats::as.formula(formula)
  vars    <- all.vars(formula)
  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)

  # keep only source rows finite on every model variable and coordinate
  fin <- Reduce(`&`, lapply(vars, function(v) is.finite(as.numeric(obs[[v]]))))
  fin <- fin & is.finite(obs$Xo) & is.finite(obs$Yo)
  obs <- obs[fin]
  if (nrow(obs) < (length(vars) + 1L)) return(NULL)

  obs_xy <- as.matrix(obs[, c("Xo", "Yo")])
  tgt_xy <- as.matrix(tgt[, c("Xt", "Yt")])
  obs_df <- as.data.frame(obs)

  if (is.null(bw)) {
    if (is.null(bw_response)) bw_response <- all.vars(formula)[1L]
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(
      coords = obs_xy, data = data.frame(.y = as.numeric(obs[[bw_response]])))
    bw <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                          kernel = kernel, adaptive = adaptive,
                          p = p, theta = theta, longlat = longlat, dMat = d_self)
    rm(d_self)
  }

  D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy,
                        focus = 0, p = p, theta = theta, longlat = longlat)
  W <- GWmodel::gw.weight(vdist = D, bw = bw, kernel = kernel, adaptive = adaptive)
  W <- matrix(W, nrow = nrow(obs_xy), ncol = nrow(tgt_xy))

  fit_one <- function(j) {
    tryCatch({
      df <- obs_df
      df[[".w"]] <- W[, j]                     # weight as a data column so lm()
      cf <- as.data.frame(stats::coef(summary( # resolves it in `data`, not env
        stats::lm(formula, data = df, weights = .w))))
      names(cf) <- c("est", "se", "tv", "pv")
      cf$term <- rownames(cf)
      cf$unit_id <- as.character(tgt[[unit]][j])
      data.table::as.data.table(cf)
    }, error = function(e) NULL)
  }

  out <- data.table::rbindlist(lapply(seq_len(nrow(tgt)), fit_one), fill = TRUE)
  if (nrow(out) == 0L) return(NULL)
  if (!is.null(terms)) out <- out[term %in% terms]
  out[, unit_level := unit]
  out[, model_estimator := "gwr_wls"]
  data.table::setcolorder(out, c("term", "unit_id", "unit_level",
                                 "est", "se", "tv", "pv", "model_estimator"))
  attr(out, "bandwidth") <- bw
  out[]
}


#' Geographically weighted local regression coefficients at points
#'
#' For each target unit, fits `lm(formula, weights = w_i)` where `w_i` are
#' `GWmodel::gw.weight()` weights over the source units, and returns the full
#' per-unit coefficient table. Accepts an arbitrary (multi-covariate) formula and
#' returns exactly what `stats::lm()` reports.
#'
#' @param data Source units (`data.table`/data frame) with the model variables
#'   and the coordinate columns in `coords`.
#' @param unit Character; the unit-id column (present in source and target).
#' @param formula A model formula, e.g. `y ~ x` or `y ~ x1 + x2`.
#' @param coords Length-2 character vector naming longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param predict_data Optional TARGET units with `unit` and `coords`. Defaults
#'   to `data` (targets = sources).
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `TRUE`.
#' @param bw Optional pre-computed bandwidth. If `NULL`, selected once via
#'   `GWmodel::bw.gwr(bw_response ~ 1)` on the source set.
#' @param bw_response Column for bandwidth selection when `bw` is `NULL`.
#'   Default: the response of `formula`.
#' @param bw_approach Bandwidth criterion. Default `"CV"`.
#' @param terms Optional character vector restricting the returned terms
#'   (e.g. `c("avail00")`). Default: all terms including the intercept.
#'
#' @return A `data.table` with one row per (target unit x term): `term`,
#'   `unit_id`, `unit_level`, `est`, `se`, `tv`, `pv`, `model_estimator`. The
#'   selected bandwidth is attached as `attr(., "bandwidth")`.
#' @family Geographically weighted regression
#' @export
estimate_gwr_coefficients_by_point <- function(data, unit, formula,
                                               coords = c("longitude", "latitude"),
                                               predict_data = NULL,
                                               distance_metric = "Euclidean",
                                               kernel = "bisquare", adaptive = TRUE,
                                               bw = NULL, bw_response = NULL,
                                               bw_approach = "CV", terms = NULL) {
  dm <- resolve_distance_metric(distance_metric)
  d  <- data.table::as.data.table(data)
  d[[unit]] <- as.character(d[[unit]])
  obs <- data.table::data.table(d)
  obs[, `:=`(Xo = as.numeric(d[[coords[1L]]]), Yo = as.numeric(d[[coords[2L]]]))]

  if (is.null(predict_data)) {
    tgt <- unique(data.table::data.table(d), by = unit)
    tgt[, `:=`(Xt = as.numeric(tgt[[coords[1L]]]), Yt = as.numeric(tgt[[coords[2L]]]))]
  } else {
    pd <- data.table::as.data.table(predict_data); pd[[unit]] <- as.character(pd[[unit]])
    tgt <- data.table::data.table(pd)
    tgt[, `:=`(Xt = as.numeric(pd[[coords[1L]]]), Yt = as.numeric(pd[[coords[2L]]]))]
  }

  .estimate_gwrcoef_core(obs = obs, tgt = tgt, unit = unit, formula = formula,
                         p = dm$p, theta = dm$theta, longlat = dm$longlat,
                         kernel = kernel, adaptive = adaptive, bw = bw,
                         bw_response = bw_response, bw_approach = bw_approach,
                         terms = terms)
}


#' Geographically weighted local regression coefficients on polygons
#'
#' As `estimate_gwr_coefficients_by_point()`, but the spatial units are polygons;
#' centroids are derived from `polygons` and joined by the unit id.
#'
#' @param data Source units with the model variables.
#' @param unit Character; the unit-id column.
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit`. Provides target ids and centroids.
#' @param formula A model formula, e.g. `y ~ x1 + x2`.
#' @param poly_id Name of the id field in `polygons`. Default: value of `unit`.
#' @param predict_ids Optional character vector of target unit ids. Default: all
#'   polygon units.
#' @param distance_metric,kernel,adaptive,bw,bw_response,bw_approach,terms
#'   See `estimate_gwr_coefficients_by_point()`.
#'
#' @return A `data.table` with one row per (target unit x term) (see
#'   `estimate_gwr_coefficients_by_point()`).
#' @family Geographically weighted regression
#' @export
estimate_gwr_coefficients_by_polygon <- function(data, unit, polygons, formula,
                                                 poly_id = unit, predict_ids = NULL,
                                                 distance_metric = "Euclidean",
                                                 kernel = "bisquare", adaptive = TRUE,
                                                 bw = NULL, bw_response = NULL,
                                                 bw_approach = "CV", terms = NULL) {
  dm  <- resolve_distance_metric(distance_metric)
  cds <- .gwkit_centroids(polygons, poly_id = poly_id, longlat = dm$longlat)

  d <- data.table::as.data.table(data)
  d[[unit]] <- as.character(d[[unit]])
  obs <- merge(d, cds, by.x = unit, by.y = "uid", all.x = FALSE)
  data.table::setnames(obs, c("X", "Y"), c("Xo", "Yo"))

  tgt <- data.table::copy(cds)
  data.table::setnames(tgt, c("uid", "X", "Y"), c(unit, "Xt", "Yt"))
  if (!is.null(predict_ids)) tgt <- tgt[get(unit) %in% as.character(predict_ids)]

  .estimate_gwrcoef_core(obs = obs, tgt = tgt, unit = unit, formula = formula,
                         p = dm$p, theta = dm$theta, longlat = dm$longlat,
                         kernel = kernel, adaptive = adaptive, bw = bw,
                         bw_response = bw_response, bw_approach = bw_approach,
                         terms = terms)
}
