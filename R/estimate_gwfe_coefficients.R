# ============================================================
# Geographically weighted FIXED-EFFECTS local regression
# ============================================================
# estimate_gwfe_coefficients_by_point()   - units carry point coordinates.
# estimate_gwfe_coefficients_by_polygon() - units are polygons; centroids derived.
#
# Localizes a panel fixed-effects model. The panel fixed effects are first
# absorbed with fixed_effect_model_data_prep() (the within transform, re-centred
# at the grand mean), so by Frisch-Waugh-Lovell the demeaned model is an ordinary
# linear one. That demeaned model is then estimated by geographically weighted
# WLS: for every spatial unit i, all panel-time observations are weighted by a
# GWmodel::gw.weight() spatial kernel over unit centroids and a weighted lm() is
# fitted. The result is a LOCAL (unit-specific) version of the national FE
# coefficients, plus a local within-R^2 that supports geographically weighted
# model selection (e.g. per-county degree-day thresholds).
#
# Returned per (unit x term): est, se, tv, pv (with an approximate FE degrees-of-
# freedom correction to the SE), and the unit-level fit columns r_squared,
# adj_r_squared, n_obs, n_units.
# ============================================================

utils::globalVariables(c("Xo", "Yo", "Xt", "Yt", "uid", "X", "Y", ".w",
                         "term", "unit_level", "model_estimator",
                         "r_squared", "adj_r_squared", "n_obs", "n_units"))

# ------------------------------------------------------------
# Internal engine. `obs` carries source coords (Xo, Yo), the panel id, and the
# (already demeaned) model variables incl. the demeaned response `dm_response`.
# `tgt` carries target coords (Xt, Yt) and target ids.
# ------------------------------------------------------------
.estimate_gwfe_core <- function(obs, tgt, unit, panel, dm_response, covs,
                                p, theta, longlat, kernel, adaptive, bw,
                                bw_approach, terms, fe_df_correction) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwfe_coefficients_*().")

  obs <- data.table::as.data.table(obs)
  tgt <- data.table::as.data.table(tgt)
  f   <- stats::as.formula(paste(dm_response, "~", paste(covs, collapse = " + ")))

  obs_xy <- as.matrix(obs[, c("Xo", "Yo")])
  tgt_xy <- as.matrix(tgt[, c("Xt", "Yt")])
  obs_df <- as.data.frame(obs)
  pcty   <- as.character(obs[[panel]])

  if (is.null(bw)) {
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(
      coords = obs_xy, data = data.frame(.y = as.numeric(obs[[dm_response]])))
    bw <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                          kernel = kernel, adaptive = adaptive,
                          p = p, theta = theta, longlat = longlat, dMat = d_self)
    rm(d_self)
  }

  D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy,
                        focus = 0, p = p, theta = theta, longlat = longlat)
  W <- GWmodel::gw.weight(vdist = D, bw = bw, kernel = kernel, adaptive = adaptive)
  W <- matrix(W, nrow = nrow(obs_xy), ncol = nrow(tgt_xy))
  ncov <- length(covs)

  fit_one <- function(j) {
    tryCatch({
      w  <- W[, j]
      df <- obs_df; df[[".w"]] <- w
      m  <- stats::lm(f, data = df, weights = .w)
      sm <- summary(m)
      cf <- as.data.frame(sm$coef)
      names(cf) <- c("est", "se", "tv", "pv")
      cf$term <- rownames(cf)

      # approximate FE df correction: rescale SE from the lm residual df to the
      # within df that also discounts the locally-absorbed unit effects.
      if (isTRUE(fe_df_correction)) {
        n_pos <- sum(w > 0)
        n_fe  <- length(unique(pcty[w > 0]))
        df_lm <- sm$df[2L]
        df_fe <- n_pos - n_fe - ncov
        if (is.finite(df_fe) && df_fe > 0 && is.finite(df_lm) && df_lm > 0) {
          cf$se <- cf$se * sqrt(df_lm / df_fe)
          cf$tv <- cf$est / cf$se
          cf$pv <- 2 * stats::pt(-abs(cf$tv), df = df_fe)
        }
      }

      cf$unit_id       <- as.character(tgt[[unit]][j])
      cf$r_squared     <- sm$r.squared
      cf$adj_r_squared <- sm$adj.r.squared
      cf$n_obs         <- sum(w > 0)
      cf$n_units       <- length(unique(pcty[w > 0]))
      data.table::as.data.table(cf)
    }, error = function(e) NULL)
  }

  out <- data.table::rbindlist(lapply(seq_len(nrow(tgt)), fit_one), fill = TRUE)
  if (nrow(out) == 0L) return(NULL)
  if (!is.null(terms)) out <- out[term %in% terms]
  out[, unit_level := unit]
  out[, model_estimator := "gwfe"]
  data.table::setcolorder(out, c("term", "unit_id", "unit_level", "est", "se",
                                 "tv", "pv", "r_squared", "adj_r_squared",
                                 "n_obs", "n_units", "model_estimator"))
  attr(out, "bandwidth") <- bw
  out[]
}

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


#' Geographically weighted fixed-effects local regression at points
#'
#' Absorbs panel fixed effects (via [fixed_effect_model_data_prep()]) and then
#' estimates the demeaned model by geographically weighted WLS, giving a local
#' (unit-specific) version of the national fixed-effects coefficients together
#' with a local within-R^2 for geographically weighted model selection.
#'
#' @param data A `data.table`/data frame panel with the model variables, the
#'   `panel` id, the `time` id, `unit`, and the coordinate columns in `coords`.
#' @param unit Character; the spatial-unit id column (present in source and
#'   target). Usually equal to `panel`.
#' @param formula A model formula, e.g. `lny ~ DD1 + DD2 + DD3 + ppt + ppt2`. The
#'   covariates are demeaned; the response is demeaned via its returned means.
#' @param panel Character vector; the fixed-effect panel id column(s).
#' @param time Character; the time id column.
#' @param coords Length-2 character vector naming longitude/latitude columns.
#'   Default `c("longitude", "latitude")`.
#' @param predict_data Optional TARGET units (one row per unit) with `unit` and
#'   `coords`. Defaults to the unique units in `data`.
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"bisquare"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `TRUE`.
#' @param bw Optional pre-computed bandwidth. If `NULL`, selected once via
#'   `GWmodel::bw.gwr()` on the demeaned response.
#' @param bw_approach Bandwidth criterion. Default `"CV"`.
#' @param terms Optional character vector restricting the returned terms.
#' @param fe_df_correction Logical; rescale the local SE from the lm residual df
#'   to a within df that discounts locally-absorbed unit effects. Default `TRUE`.
#'
#' @return A `data.table` with one row per (unit x term): `term`, `unit_id`,
#'   `unit_level`, `est`, `se`, `tv`, `pv`, `r_squared`, `adj_r_squared`,
#'   `n_obs`, `n_units`, `model_estimator` (`"gwfe"`). Bandwidth in
#'   `attr(., "bandwidth")`; number of panels absorbed in `attr(., "NFE")`.
#' @family Geographically weighted regression
#' @export
estimate_gwfe_coefficients_by_point <- function(data, unit, formula, panel, time,
                                                coords = c("longitude", "latitude"),
                                                predict_data = NULL,
                                                distance_metric = "Euclidean",
                                                kernel = "bisquare", adaptive = TRUE,
                                                bw = NULL, bw_approach = "CV",
                                                terms = NULL, fe_df_correction = TRUE) {
  dm  <- resolve_distance_metric(distance_metric)
  dem <- .gwfe_demean(data, formula, panel, time)
  d   <- dem$data
  d[[unit]] <- as.character(d[[unit]])

  # the demean step keeps only panel/time/model columns, so re-attach the
  # unit-level (constant-within-unit) coordinates from the original data.
  src <- data.table::as.data.table(data); src[[unit]] <- as.character(src[[unit]])
  coord_map <- unique(src[, c(unit, coords), with = FALSE], by = unit)
  d <- merge(d, coord_map, by = unit, all.x = TRUE)

  obs <- data.table::copy(d)
  obs[, `:=`(Xo = as.numeric(d[[coords[1L]]]), Yo = as.numeric(d[[coords[2L]]]))]

  if (is.null(predict_data)) {
    tgt <- unique(data.table::copy(d), by = unit)
    tgt[, `:=`(Xt = as.numeric(tgt[[coords[1L]]]), Yt = as.numeric(tgt[[coords[2L]]]))]
  } else {
    pd <- data.table::as.data.table(predict_data); pd[[unit]] <- as.character(pd[[unit]])
    tgt <- pd; tgt[, `:=`(Xt = as.numeric(pd[[coords[1L]]]), Yt = as.numeric(pd[[coords[2L]]]))]
  }

  out <- .estimate_gwfe_core(obs, tgt, unit = unit, panel = panel,
                             dm_response = dem$dm_response, covs = dem$covs,
                             p = dm$p, theta = dm$theta, longlat = dm$longlat,
                             kernel = kernel, adaptive = adaptive, bw = bw,
                             bw_approach = bw_approach, terms = terms,
                             fe_df_correction = fe_df_correction)
  if (!is.null(out)) attr(out, "NFE") <- dem$NFE
  out
}


#' Geographically weighted fixed-effects local regression on polygons
#'
#' As `estimate_gwfe_coefficients_by_point()`, but the spatial units are
#' polygons; centroids are derived from `polygons` and joined by the unit id.
#'
#' @param data A panel `data.table`/data frame with the model variables, the
#'   `panel` id and `time` id.
#' @param unit Character; the spatial-unit id column.
#' @param polygons A `terra` `SpatVector` or `sf` polygon layer whose id field
#'   (`poly_id`) matches `unit`.
#' @param formula,panel,time See `estimate_gwfe_coefficients_by_point()`.
#' @param poly_id Name of the id field in `polygons`. Default: value of `unit`.
#' @param predict_ids Optional character vector of target unit ids.
#' @param distance_metric,kernel,adaptive,bw,bw_approach,terms,fe_df_correction
#'   See `estimate_gwfe_coefficients_by_point()`.
#'
#' @return A `data.table` (see `estimate_gwfe_coefficients_by_point()`).
#' @family Geographically weighted regression
#' @export
estimate_gwfe_coefficients_by_polygon <- function(data, unit, polygons, formula,
                                                  panel, time, poly_id = unit,
                                                  predict_ids = NULL,
                                                  distance_metric = "Euclidean",
                                                  kernel = "bisquare", adaptive = TRUE,
                                                  bw = NULL, bw_approach = "CV",
                                                  terms = NULL, fe_df_correction = TRUE) {
  dm  <- resolve_distance_metric(distance_metric)
  dem <- .gwfe_demean(data, formula, panel, time)
  d   <- dem$data
  d[[unit]] <- as.character(d[[unit]])

  cds <- .gwkit_centroids(polygons, poly_id = poly_id, longlat = dm$longlat)
  obs <- merge(d, cds, by.x = unit, by.y = "uid", all.x = FALSE)
  data.table::setnames(obs, c("X", "Y"), c("Xo", "Yo"))

  tgt <- data.table::copy(cds)
  data.table::setnames(tgt, c("uid", "X", "Y"), c(unit, "Xt", "Yt"))
  if (!is.null(predict_ids)) tgt <- tgt[get(unit) %in% as.character(predict_ids)]

  out <- .estimate_gwfe_core(obs, tgt, unit = unit, panel = panel,
                             dm_response = dem$dm_response, covs = dem$covs,
                             p = dm$p, theta = dm$theta, longlat = dm$longlat,
                             kernel = kernel, adaptive = adaptive, bw = bw,
                             bw_approach = bw_approach, terms = terms,
                             fe_df_correction = fe_df_correction)
  if (!is.null(out)) attr(out, "NFE") <- dem$NFE
  out
}
