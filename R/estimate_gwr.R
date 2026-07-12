# ============================================================
# Geographically weighted local regression (unified)
# ============================================================
# estimate_gwr() - one entry point; point vs polygon is detected by class (an
# sf/SpatVector input -> polygon units reduced to point-on-surface; otherwise
# `coords` columns are used) via .resolve_gw_geometry().
#
# ONE geographically weighted regression, mode-switched by arguments, all over
# the shared .gw_local_fit() engine (GWmodel::gw.dist + GWmodel::gw.weight +
# stats::lm.wfit). It subsumes the former estimate_gwr / estimate_gwr_coefficients
# / estimate_gwfe_coefficients:
#
#   * formula          - arbitrary (one or many covariates).
#   * panel + time      - if given, panel FIXED EFFECTS are absorbed first via
#                         fixed_effect_model_data_prep() (Frisch-Waugh-Lovell),
#                         so the demeaned model is a plain WLS. model_estimator
#                         becomes "gwfe"; otherwise "gwr".
#   * variance = TRUE   - also fit the local slope of log-squared residuals (a
#                         heteroskedasticity companion), returned as extra rows.
#   * fit_stats = TRUE  - attach the local within-R^2 and n_obs (model/knot
#                         selection).
#   * standard_errors   - analytic WLS SEs (default FALSE; inference is meant to
#                         come from the bootstrap).
#
# Output is LONG: one row per (unit_id x term x estimand), estimand in
# {"mean","variance"}, columns estimate/se/tv/pv (+ r_squared/n_obs when
# fit_stats). Point estimates are identical to a per-unit stats::lm().
# ============================================================

utils::globalVariables(c("Xo", "Yo", "Xt", "Yt", "uid", "X", "Y",
                         "longitude", "latitude", "term", "estimand", "estimate",
                         "se", "tv", "pv", "r_squared", "n_obs",
                         "unit_level", "model_estimator"))


# ------------------------------------------------------------
# Internal: demean a panel with fixed_effect_model_data_prep() and return the
# demeaned data.table (with a demeaned response column) plus the response name.
# ------------------------------------------------------------
.gwfe_demean <- function(data, formula, panel, time) {
  formula  <- stats::as.formula(formula)
  response <- all.vars(formula)[1L]
  covs     <- attr(stats::terms(formula), "term.labels")
  prep <- fixed_effect_model_data_prep(data = data, varlist = covs, panel = panel,
                                       time = time, output = response)
  d <- prep$data
  dm_response <- paste0(response, "_dm")
  d[[dm_response]] <- d[[response]] - d[[paste0(response, "_mean_i")]] +
                                       d[[paste0(response, "_mean")]]
  list(data = d, covs = covs, dm_response = dm_response, NFE = prep$NFE)
}


# ------------------------------------------------------------
# Internal engine. `obs` carries source coords (Xo, Yo) and the (already
# demeaned, if FE) model variables; `tgt` carries target coords (Xt, Yt) and ids.
# ------------------------------------------------------------
.estimate_gwr_core <- function(obs, tgt, unit, formula,
                               p, theta, longlat, kernel, adaptive, bw, bw_approach,
                               variance, fit_stats, standard_errors, terms, eps,
                               model_estimator) {

  if (!kernel %in% c("gaussian", "exponential", "bisquare", "boxcar", "tricube"))
    stop("Unknown kernel: ", kernel)

  formula <- stats::as.formula(formula)
  vars <- all.vars(formula)
  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)

  fin <- Reduce(`&`, lapply(vars, function(v) is.finite(as.numeric(obs[[v]]))))
  fin <- fin & is.finite(obs$Xo) & is.finite(obs$Yo)
  obs <- obs[fin]
  if (nrow(obs) < (length(vars) + 1L)) return(NULL)

  obs_df <- as.data.frame(obs)
  mf <- stats::model.frame(formula, obs_df)
  mm <- stats::model.matrix(attr(mf, "terms"), mf)   # built once (incl. intercept)
  y  <- as.numeric(stats::model.response(mf))
  obs_xy <- as.matrix(obs[, c("Xo", "Yo")])
  tgt_xy <- as.matrix(tgt[, c("Xt", "Yt")])

  fit <- .gw_local_fit(mm = mm, y = y, obs_xy = obs_xy, tgt_xy = tgt_xy,
                       p = p, theta = theta, longlat = longlat,
                       kernel = kernel, adaptive = adaptive, bw = bw,
                       bw_approach = bw_approach,
                       standard_errors = standard_errors, fit_stats = fit_stats)
  bw <- fit$bw
  terms_all <- colnames(mm)
  uid <- as.character(tgt[[unit]])

  long_rows <- function(f, estimand) {
    data.table::rbindlist(lapply(seq_along(terms_all), function(ti){
      dt <- data.table::data.table(
        unit_id  = uid,
        term     = terms_all[ti],
        estimand = estimand,
        estimate = f$coef[, ti],
        se       = if (!is.null(f$se)) f$se[, ti] else NA_real_,
        tv = NA_real_, pv = NA_real_)
      if (isTRUE(fit_stats)) { dt[, r_squared := f$r_squared]; dt[, n_obs := f$n_obs] }
      dt
    }), fill = TRUE)
  }

  out <- long_rows(fit, "mean")

  if (isTRUE(variance)) {
    ridx  <- match(as.character(obs[[unit]]), uid)          # obs -> its unit's fit
    fitted_obs <- rowSums(mm * fit$coef[ridx, , drop = FALSE])
    lrs   <- log((y - fitted_obs)^2 + eps)
    vfit  <- .gw_local_fit(mm = mm, y = lrs, obs_xy = obs_xy, tgt_xy = tgt_xy,
                           p = p, theta = theta, longlat = longlat,
                           kernel = kernel, adaptive = adaptive, bw = bw,
                           bw_approach = bw_approach,
                           standard_errors = standard_errors, fit_stats = fit_stats)
    out <- data.table::rbindlist(list(out, long_rows(vfit, "variance")), fill = TRUE)
  }

  if (isTRUE(standard_errors)) {
    out[, tv := estimate / se]
    out[, pv := 2 * stats::pnorm(-abs(tv))]
  }
  if (!is.null(terms)) out <- out[term %in% terms]
  out[, unit_level := unit]
  out[, model_estimator := model_estimator]
  data.table::setcolorder(out, c("unit_id", "term", "estimand", "unit_level",
                                 "estimate", "se", "tv", "pv"))
  attr(out, "bandwidth") <- bw
  out[]
}


#' Geographically weighted local regression
#'
#' One geographically weighted regression covering plain GWR, fixed-effects GWR,
#' and the mean/variance decomposition, over the shared local WLS engine. Point
#' vs polygon mode is detected by class: pass an `sf`/`terra` `SpatVector` (as
#' `data`, or via `geometry` for panel inputs) for polygon units - each reduced
#' to its point-on-surface; otherwise the `coords` columns of `data` are used.
#'
#' @param data A `data.table`/data frame with the model variables (and, for
#'   fixed effects, the `panel`/`time` columns), OR an `sf`/`SpatVector` carrying
#'   both attributes and geometry (polygon mode).
#' @param unit Character; the spatial-unit id column (equal to `panel` in
#'   fixed-effects mode).
#' @param formula A model formula, e.g. `y ~ x` or `y ~ x1 + x2`.
#' @param geometry Optional `sf`/`SpatVector` polygon layer (whose `poly_id`
#'   matches `unit`) supplying per-unit geometry when `data` is a plain (e.g.
#'   panel) table. If `NULL` and `data` is not spatial, point mode is used.
#' @param coords Length-2 character vector naming longitude/latitude columns
#'   (point mode). Default `c("longitude", "latitude")`.
#' @param poly_id Name of the id field in the polygon layer. Default: `unit`.
#' @param predict Optional target restriction: a data frame of point targets
#'   (point mode) or a character id vector / data frame with `unit` (polygon
#'   mode). Defaults to the unique units in `data`.
#' @param panel,time Optional fixed-effects panel id and time id. If `panel` is
#'   supplied the panel fixed effects are absorbed (within transform) before the
#'   GW fit and `model_estimator` is `"gwfe"`; otherwise `"gwr"`.
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `TRUE`.
#' @param bw Optional pre-computed bandwidth. If `NULL`, selected once via
#'   `GWmodel::bw.gwr()` on the (demeaned) response.
#' @param bw_approach Bandwidth criterion. Default `"CV"`.
#' @param variance Logical; also fit the local slope of the log-squared residuals
#'   (a heteroskedasticity companion), returned as `estimand == "variance"` rows.
#'   Requires the targets to cover the source units. Default `FALSE`.
#' @param fit_stats Logical; attach local within-`r_squared` and `n_obs`.
#'   Default `FALSE`.
#' @param standard_errors Logical; also return analytic WLS `se`/`tv`/`pv`.
#'   Default `FALSE` - coefficients only (identical to `lm`); use the bootstrap.
#' @param terms Optional character vector restricting the returned terms.
#' @param eps Added to squared residuals before the log (variance mode). Default
#'   `1e-8`.
#'
#' @return A `data.table`, one row per `unit_id` x `term` x `estimand`
#'   (`estimand` in `"mean"`/`"variance"`): `unit_id`, `term`, `estimand`,
#'   `unit_level`, `estimate`, `se`, `tv`, `pv`, (when `fit_stats`) `r_squared`,
#'   `n_obs`, and `model_estimator`. Bandwidth in `attr(., "bandwidth")`; number
#'   of panels absorbed (FE mode) in `attr(., "NFE")`.
#' @family Geographically weighted regression
#' @export
estimate_gwr <- function(data, unit, formula, geometry = NULL,
                         coords = c("longitude", "latitude"), poly_id = unit,
                         predict = NULL, panel = NULL, time = NULL,
                         distance_metric = "Euclidean", kernel = "bisquare",
                         adaptive = TRUE, bw = NULL, bw_approach = "CV",
                         variance = FALSE, fit_stats = FALSE,
                         standard_errors = FALSE, terms = NULL, eps = 1e-8) {
  dm <- resolve_distance_metric(distance_metric)
  formula <- stats::as.formula(formula)
  nfe <- NA_integer_
  is_poly <- inherits(geometry, c("sf", "SpatVector")) ||
             inherits(data, c("sf", "SpatVector"))

  if (!is.null(panel)) {
    dem <- .gwfe_demean(data, formula, panel, time)
    d   <- dem$data; d[[unit]] <- as.character(d[[unit]])
    formula <- stats::reformulate(dem$covs, response = dem$dm_response)
    model_est <- "gwfe"; nfe <- dem$NFE
    if (!is_poly) {                       # point + FE: re-attach coords by unit
      src  <- data.table::as.data.table(data); src[[unit]] <- as.character(src[[unit]])
      cmap <- unique(src[, c(unit, coords), with = FALSE], by = unit)
      d <- merge(d, cmap, by = unit, all.x = TRUE)
    }
    data_for_geo <- d
  } else {
    model_est <- "gwr"; data_for_geo <- data
  }

  geo <- .resolve_gw_geometry(data_for_geo, unit = unit, coords = coords,
                              geometry = geometry, poly_id = poly_id,
                              longlat = dm$longlat, predict = predict)

  out <- .estimate_gwr_core(geo$obs, geo$tgt, unit = unit, formula = formula,
                            p = dm$p, theta = dm$theta, longlat = dm$longlat,
                            kernel = kernel, adaptive = adaptive, bw = bw,
                            bw_approach = bw_approach, variance = variance,
                            fit_stats = fit_stats, standard_errors = standard_errors,
                            terms = terms, eps = eps, model_estimator = model_est)
  if (!is.null(out)) attr(out, "NFE") <- nfe
  out
}
