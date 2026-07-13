# ============================================================
# Geographically weighted summary statistics (GWSS)
# ============================================================
# estimate_gwss() - one entry point over GWmodel::gwss(). At each prediction
# location it returns the local summary statistics (local mean, SD, CV, skew, and
# local covariance/correlation between variables). Point vs polygon mode is
# class-detected:
#
#   * point mode  - `data` is a plain table with coordinate columns (`coords`)
#                   and `variable_list`; prediction locations are `predict`
#                   (default: the observed points), radius-screened.
#   * polygon mode- `data` is joined by `id_col` to an sf/SpatVector `geometry`;
#                   each polygon is reduced to its point-on-surface, GWSS is fit
#                   on the observed (finite-value) polygons and evaluated at ALL
#                   polygons (spatial smoothing / gap-filling).
#   * us_counties = TRUE pulls the county geometry from `urbnmapr` (a convenience
#                   polygon mode; `id_col` defaults to "county_fips").
#
# Bandwidth is selected once (GWmodel::bw.gwr on a random observed subsample);
# distances/weights are GWmodel's. Coordinates are assumed WGS84 (EPSG:4326) in
# point mode; polygons are projected to `target_crs` (or 4326 for great-circle).
# ============================================================

utils::globalVariables(c("value", "longitude", "latitude",
                         "longitude_flag", "latitude_flag", ".uid"))


# ------------------------------------------------------------
# Internal: metric-dependent prep. Builds the seeded bandwidth subsample and the
# two distance matrices (obs-obs subsample + obs-predict) that do NOT depend on
# the kernel, so they can be reused across kernels for the same distance metric.
# `pts_sp_obs`/`pts_sp_all` are SpatialPointsDataFrames; `pts_sp_obs@data` must
# carry a numeric `value` column. Returns list(dMat_sub, dMat_os, sub_ids, n_obs)
# or NULL. Callers must snapshot/restore the RNG (this calls set.seed).
# ------------------------------------------------------------
.gwss_prep <- function(pts_sp_obs, pts_sp_all, p, theta, longlat, draw_rate) {
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required for estimate_gwss().")

  coords_obs <- sp::coordinates(pts_sp_obs)
  coords_all <- sp::coordinates(pts_sp_all)
  n_obs <- nrow(pts_sp_obs)
  if (n_obs < 5L) { message("Too few observed units to fit GW smoothing."); return(NULL) }

  n_sub <- min(n_obs - 1L, max(5L, ceiling(draw_rate * n_obs)))
  set.seed(1L)
  sub_ids    <- sample.int(n_obs, n_sub)
  coords_sub <- coords_obs[sub_ids, , drop = FALSE]

  dMat_sub <- GWmodel::gw.dist(dp.locat = coords_sub, rp.locat = coords_sub,
                               p = p, theta = theta, longlat = longlat)
  dMat_os  <- GWmodel::gw.dist(dp.locat = coords_obs, rp.locat = coords_all,
                               p = p, theta = theta, longlat = longlat)
  list(dMat_sub = dMat_sub, dMat_os = dMat_os, sub_ids = sub_ids, n_obs = n_obs)
}

# ------------------------------------------------------------
# Internal: kernel-dependent run, reusing the cached `prep` matrices. If `bw` is
# NULL the bandwidth is selected (bw.gwr on the subsample) for THIS kernel;
# pass a numeric `bw` to reuse a bandwidth across kernels/metrics. Returns
# list(df, bw) or NULL.
# ------------------------------------------------------------
.gwss_run <- function(pts_sp_obs, pts_sp_all, variable_list, prep,
                      kernel, adaptive, approach, p, theta, longlat, bw = NULL) {
  if (is.null(bw))
    bw <- GWmodel::bw.gwr(value ~ 1, data = pts_sp_obs[prep$sub_ids, ], approach = approach,
                          adaptive = adaptive, kernel = kernel,
                          p = p, theta = theta, longlat = longlat, dMat = prep$dMat_sub)
  gwss_obj <- GWmodel::gwss(data = pts_sp_obs, summary.locat = pts_sp_all,
                            bw = bw, vars = variable_list, kernel = kernel,
                            adaptive = adaptive, p = p, theta = theta,
                            longlat = longlat, dMat = prep$dMat_os, quantile = FALSE)
  list(df = data.table::as.data.table(as.data.frame(gwss_obj$SDF@data)), bw = bw)
}

# ------------------------------------------------------------
# Internal: original single-spec fit (prep + run), RNG-guarded. Kept so the
# scalar estimate_gwss() path is byte-for-byte unchanged.
# ------------------------------------------------------------
.gwss_fit <- function(pts_sp_obs, pts_sp_all, variable_list,
                      p, theta, longlat, kernel, adaptive, approach, draw_rate) {
  # Deterministic subsample WITHOUT disturbing the caller's RNG stream.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  prep <- .gwss_prep(pts_sp_obs, pts_sp_all, p, theta, longlat, draw_rate)
  if (is.null(prep)) return(NULL)
  .gwss_run(pts_sp_obs, pts_sp_all, variable_list, prep,
            kernel, adaptive, approach, p, theta, longlat)
}

# ------------------------------------------------------------
# Internal: point-mode setup shared by estimate_gwss() and
# estimate_gwss_kernels(). Builds the observed/prediction SpatialPointsDataFrames
# and applies the metric-independent radius pre-screen. Returns
# list(obs_sp, pred_sp); `pred_sp` may have 0 rows (nothing passed the screen).
# ------------------------------------------------------------
.gwss_point_setup <- function(data, variable_list, coords, predict,
                              target_crs, feasible_radius) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

  obs <- data.table::copy(data)
  data.table::setnames(obs, coords, c("longitude", "latitude"))
  obs[, value := as.numeric(get(variable_list[1]))]
  obs <- obs[is.finite(value) & is.finite(longitude) & is.finite(latitude)]
  obs <- obs[, lapply(.SD, mean, na.rm = TRUE), by = .(longitude, latitude),
             .SDcols = unique(c("value", variable_list))]

  if (is.null(predict)) {
    pred <- data.table::copy(obs)
  } else {
    pred <- data.table::as.data.table(predict)
    data.table::setnames(pred, coords, c("longitude", "latitude"), skip_absent = TRUE)
  }
  pred <- pred[is.finite(longitude) & is.finite(latitude)]
  pred[, longitude_flag := longitude]; pred[, latitude_flag := latitude]

  obs_sp  <- as.data.frame(obs);  sp::coordinates(obs_sp)  <- ~ longitude + latitude
  pred_sp <- as.data.frame(pred); sp::coordinates(pred_sp) <- ~ longitude + latitude

  # radius pre-screen (projected metric CRS) - independent of metric/kernel
  obs_sf <- sf::st_as_sf(obs_sp);  sf::st_crs(obs_sf) <- 4326
  all_sf <- sf::st_as_sf(pred_sp); sf::st_crs(all_sf) <- 4326
  idx <- sf::st_is_within_distance(sf::st_transform(all_sf, target_crs),
                                   sf::st_transform(obs_sf, target_crs),
                                   dist = feasible_radius * 1609.344)
  pred_sp <- pred_sp[lengths(idx) > 0, ]
  list(obs_sp = obs_sp, pred_sp = pred_sp)
}

# ------------------------------------------------------------
# Internal: attach point-mode coordinates + identifiers to a GWSS result frame.
# ------------------------------------------------------------
.gwss_point_attach <- function(gw_df, pred_sp, identifiers) {
  gw_df[, longitude := pred_sp@data[["longitude_flag"]]]
  gw_df[, latitude  := pred_sp@data[["latitude_flag"]]]
  if (!is.null(identifiers))
    for (ix in identifiers) gw_df[[ix]] <- pred_sp@data[[ix]]
  gw_df
}


#' Estimate geographically weighted summary statistics (GWSS)
#'
#' One entry point over `GWmodel::gwss()`: at each prediction location it returns
#' local summary statistics (local mean, SD, CV, skewness, and local covariance /
#' correlation between variables). Point vs polygon mode is detected by class.
#'
#' @param data A `data.table`/data frame with `variable_list` and either the
#'   `coords` columns (point mode) or the `id_col` join key (polygon mode).
#' @param variable_list Character vector of numeric variable column(s) to
#'   summarize (passed to `GWmodel::gwss(vars = )`). The bandwidth uses the first.
#' @param geometry Optional `sf`/`SpatVector` polygon layer whose `id_col` field
#'   matches `data[[id_col]]` (polygon mode). If `NULL` and `us_counties = FALSE`,
#'   point mode is used.
#' @param us_counties Logical; if `TRUE`, pull the U.S. county polygons from
#'   `urbnmapr` (polygon mode). `id_col` defaults to `"county_fips"`.
#' @param id_col Character; the unit id column shared by `data` and `geometry`
#'   (polygon mode). Required in polygon mode.
#' @param coords Length-2 character vector naming longitude/latitude columns
#'   (point mode). Default `c("longitude", "latitude")`.
#' @param predict Optional prediction locations (point mode): a data frame with
#'   the `coords` columns (and any `identifiers`). Defaults to the observed points.
#' @param distance_metric One of `gw_distance_metric_names()`. Default
#'   `"Euclidean"`.
#' @param kernel GW kernel. Default `"gaussian"`.
#' @param adaptive Logical; adaptive (kNN) bandwidth if `TRUE`. Default `TRUE`.
#' @param approach Bandwidth criterion, one of `"CV"`, `"AIC"`, `"AICc"`.
#'   Default `"CV"`.
#' @param target_crs Integer EPSG code (meters) for projection (polygon geometry
#'   and the point-mode radius screen). Default `5070` (NAD83 / CONUS Albers).
#' @param draw_rate Fraction of observed units used for bandwidth selection.
#'   Default `0.5`.
#' @param feasible_radius Radius in miles for pre-screening prediction points
#'   (point mode only). Default `100`.
#' @param identifiers Optional character vector of columns in `predict` (point
#'   mode) to carry through to the output.
#' @param state_fips_limits Optional state FIPS to restrict the county set
#'   (`us_counties`/polygon mode, when a `state_fips` column is present).
#'
#' @return A `data.table` of GWSS statistics with one row per prediction unit,
#'   the unit id (polygon mode) or coordinates + `identifiers` (point mode). The
#'   bandwidth and distance parameters are attached as attributes.
#' @family Geographically weighted summaries
#' @export
estimate_gwss <- function(data, variable_list, geometry = NULL, us_counties = FALSE,
                          id_col = NULL, coords = c("longitude", "latitude"),
                          predict = NULL, distance_metric = "Euclidean",
                          kernel = "gaussian", adaptive = TRUE, approach = "CV",
                          target_crs = 5070, draw_rate = 0.5, feasible_radius = 100,
                          identifiers = NULL, state_fips_limits = NULL) {

  if (!kernel %in% c("gaussian", "exponential", "bisquare", "boxcar", "tricube"))
    stop("`kernel` must be one of: gaussian, exponential, bisquare, boxcar, tricube")
  if (!approach %in% c("CV", "AIC", "AICc"))
    stop("`approach` must be one of: CV, AIC, AICc")
  if (missing(data) || is.null(data)) stop("Argument `data` must be supplied.")

  dm <- resolve_distance_metric(distance_metric)
  p <- dm$p; theta <- dm$theta; longlat <- dm$longlat
  data <- data.table::as.data.table(data)

  if (isTRUE(us_counties)) {
    if (!requireNamespace("urbnmapr", quietly = TRUE))
      stop("Package 'urbnmapr' is required for `us_counties = TRUE`.")
    geometry <- urbnmapr::get_urbn_map("counties", sf = TRUE)
    if (is.null(id_col)) id_col <- "county_fips"
  }

  is_poly <- inherits(geometry, c("sf", "SpatVector"))

  # ---- Polygon mode (supplied geometry or us_counties) -----------------------
  if (is_poly) {
    if (is.null(id_col)) stop("`id_col` is required in polygon mode.")
    if (!all(c(id_col, variable_list) %in% names(data)))
      stop("`data` must contain: ", id_col, " and ", paste(variable_list, collapse = ", "))
    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

    gsf <- if (inherits(geometry, "SpatVector")) sf::st_as_sf(geometry) else geometry
    if (!id_col %in% names(gsf)) stop("`geometry` must contain the `", id_col, "` column.")
    gsf$.uid <- as.character(gsf[[id_col]])
    gsf <- if (isTRUE(longlat)) sf::st_transform(gsf, 4326) else sf::st_transform(gsf, target_crs)

    d <- data.table::copy(data)
    d[, .uid := as.character(get(id_col))]
    d[, value := as.numeric(get(variable_list[1]))]
    d <- d[is.finite(value) & !is.na(.uid) & nzchar(.uid)]
    d <- unique(d, by = ".uid")

    sf_join <- dplyr::left_join(gsf, as.data.frame(d), by = ".uid")
    if (!is.null(state_fips_limits) && "state_fips" %in% names(sf_join))
      sf_join <- sf_join[as.numeric(as.character(sf_join$state_fips)) %in%
                           as.numeric(as.character(state_fips_limits)), ]
    sf_obs <- dplyr::filter(sf_join, is.finite(value))
    if (nrow(sf_obs) < 5L) { message("Too few observed units to fit GW smoothing."); return(NULL) }

    pts_all <- methods::as(suppressWarnings(sf::st_point_on_surface(sf_join)), "Spatial")
    pts_obs <- methods::as(suppressWarnings(sf::st_point_on_surface(sf_obs)),  "Spatial")

    fit <- .gwss_fit(pts_obs, pts_all, variable_list, p, theta, longlat,
                     kernel, adaptive, approach, draw_rate)
    if (is.null(fit)) return(NULL)
    gw_df <- fit$df
    gw_df[[id_col]] <- pts_all@data[[".uid"]]

  # ---- Point mode ------------------------------------------------------------
  } else {
    if (!is.null(geometry))
      stop("`geometry` must be a terra SpatVector or an sf object.")
    if (!all(c(coords, variable_list) %in% names(data)))
      stop("`data` must contain: ", paste(c(coords, variable_list), collapse = ", "))

    setup <- .gwss_point_setup(data, variable_list, coords, predict,
                               target_crs, feasible_radius)
    if (nrow(setup$pred_sp) == 0L) return(data.table::data.table())

    fit <- .gwss_fit(setup$obs_sp, setup$pred_sp, variable_list, p, theta, longlat,
                     kernel, adaptive, approach, draw_rate)
    if (is.null(fit)) return(NULL)
    gw_df <- .gwss_point_attach(fit$df, setup$pred_sp, identifiers)
  }

  data.table::setattr(gw_df, "bandwidth", fit$bw)
  data.table::setattr(gw_df, "distance_params", list(p = p, theta = theta, longlat = longlat))
  data.table::setattr(gw_df, "kernel", kernel)
  data.table::setattr(gw_df, "approach", approach)
  data.table::setattr(gw_df, "adaptive", adaptive)
  gw_df[]
}


#' Estimate GWSS across kernels for one distance metric (point mode)
#'
#' Point-mode companion to `estimate_gwss()` that evaluates several kernels for a
#' single `distance_metric` while computing the (kernel-independent) distance
#' matrices and radius screen only **once**, reusing them across kernels. This is
#' the per-metric building block for spec ensembles: much cheaper than calling
#' `estimate_gwss()` separately for each kernel.
#'
#' @inheritParams estimate_gwss
#' @param kernel Character vector of kernels to evaluate (any of `"gaussian"`,
#'   `"exponential"`, `"bisquare"`, `"boxcar"`, `"tricube"`).
#' @param bandwidth Bandwidth strategy across the supplied kernels: `"shared"`
#'   (default) selects one bandwidth (from the first kernel) and reuses it for all
#'   kernels - one `bw.gwr()` call instead of `length(kernel)`; `"per_kernel"`
#'   selects a bandwidth per kernel (matches looping `estimate_gwss()` exactly).
#' @param bw Optional numeric bandwidth override. When supplied it is used for
#'   every kernel and no bandwidth selection is performed (overrides `bandwidth`).
#'
#' @return A `data.table` stacked over `kernel` (one block of prediction rows per
#'   kernel), with the GWSS statistics, coordinates, `identifiers`, and a `kernel`
#'   column. `distance_params` is attached as an attribute; when a single shared
#'   bandwidth is used it is attached as `bandwidth`.
#' @family Geographically weighted summaries
#' @export
estimate_gwss_kernels <- function(data, variable_list, coords = c("longitude", "latitude"),
                                  predict = NULL, distance_metric = "Euclidean",
                                  kernel = c("gaussian", "exponential", "bisquare",
                                             "boxcar", "tricube"),
                                  adaptive = TRUE, approach = "CV", target_crs = 5070,
                                  draw_rate = 0.5, feasible_radius = 100,
                                  identifiers = NULL,
                                  bandwidth = c("shared", "per_kernel"), bw = NULL) {
  bandwidth <- match.arg(bandwidth)
  ok_kernels <- c("gaussian", "exponential", "bisquare", "boxcar", "tricube")
  bad_k <- setdiff(kernel, ok_kernels)
  if (length(bad_k)) stop("Unknown kernel(s): ", paste(bad_k, collapse = ", "))
  if (!approach %in% c("CV", "AIC", "AICc"))
    stop("`approach` must be one of: CV, AIC, AICc")
  if (missing(data) || is.null(data)) stop("Argument `data` must be supplied.")

  dm <- resolve_distance_metric(distance_metric)
  p <- dm$p; theta <- dm$theta; longlat <- dm$longlat
  data <- data.table::as.data.table(data)

  setup <- .gwss_point_setup(data, variable_list, coords, predict,
                             target_crs, feasible_radius)
  if (nrow(setup$pred_sp) == 0L) return(data.table::data.table())

  # Deterministic subsample WITHOUT disturbing the caller's RNG stream.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }

  prep <- .gwss_prep(setup$obs_sp, setup$pred_sp, p, theta, longlat, draw_rate)
  if (is.null(prep)) return(NULL)

  # Bandwidth reuse: an explicit `bw` overrides everything; otherwise "shared"
  # selects it once (first kernel) and reuses it across all kernels.
  shared_bw <- bw
  if (is.null(shared_bw) && bandwidth == "shared")
    shared_bw <- GWmodel::bw.gwr(value ~ 1, data = setup$obs_sp[prep$sub_ids, ],
                                 approach = approach, adaptive = adaptive,
                                 kernel = kernel[[1]], p = p, theta = theta,
                                 longlat = longlat, dMat = prep$dMat_sub)

  parts <- lapply(kernel, function(km) {
    r <- .gwss_run(setup$obs_sp, setup$pred_sp, variable_list, prep,
                   km, adaptive, approach, p, theta, longlat, bw = shared_bw)
    if (is.null(r)) return(NULL)
    d <- .gwss_point_attach(r$df, setup$pred_sp, identifiers)
    d[, kernel := km]
    d
  })

  out <- data.table::rbindlist(Filter(Negate(is.null), parts),
                               use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(data.table::data.table())
  data.table::setattr(out, "distance_params", list(p = p, theta = theta, longlat = longlat))
  if (!is.null(shared_bw)) data.table::setattr(out, "bandwidth", shared_bw)
  out[]
}


#' Select a geographically weighted bandwidth (point mode)
#'
#' Lightweight bandwidth-only helper: runs the point-mode setup + radius screen
#' and selects a single GW bandwidth via `GWmodel::bw.gwr()` on the seeded
#' subsample, **without** computing any GWSS surfaces. Useful for reusing one
#' bandwidth across many specs (e.g. a "shared" bandwidth across metrics/kernels).
#' The value matches what `estimate_gwss()` / `estimate_gwss_kernels()` would
#' select for the same `data`, `distance_metric`, `kernel`, and `draw_rate`.
#'
#' @inheritParams estimate_gwss
#' @param kernel A single GW kernel used for bandwidth selection. Default
#'   `"gaussian"`.
#'
#' @return A single numeric bandwidth (adaptive kNN count when `adaptive = TRUE`,
#'   else a distance), or `NA_real_` if no prediction point passes the radius
#'   screen or there are too few observed units.
#' @family Geographically weighted summaries
#' @export
gw_select_bandwidth <- function(data, variable_list, coords = c("longitude", "latitude"),
                                predict = NULL, distance_metric = "Euclidean",
                                kernel = "gaussian", adaptive = TRUE, approach = "CV",
                                target_crs = 5070, draw_rate = 0.5, feasible_radius = 100) {
  ok_kernels <- c("gaussian", "exponential", "bisquare", "boxcar", "tricube")
  if (length(kernel) != 1L || !kernel %in% ok_kernels)
    stop("`kernel` must be a single one of: ", paste(ok_kernels, collapse = ", "))
  if (!approach %in% c("CV", "AIC", "AICc"))
    stop("`approach` must be one of: CV, AIC, AICc")
  if (missing(data) || is.null(data)) stop("Argument `data` must be supplied.")

  dm <- resolve_distance_metric(distance_metric)
  p <- dm$p; theta <- dm$theta; longlat <- dm$longlat
  data <- data.table::as.data.table(data)

  setup <- .gwss_point_setup(data, variable_list, coords, predict,
                             target_crs, feasible_radius)
  if (nrow(setup$pred_sp) == 0L) return(NA_real_)

  # Deterministic subsample WITHOUT disturbing the caller's RNG stream.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }

  prep <- .gwss_prep(setup$obs_sp, setup$pred_sp, p, theta, longlat, draw_rate)
  if (is.null(prep)) return(NA_real_)

  GWmodel::bw.gwr(value ~ 1, data = setup$obs_sp[prep$sub_ids, ], approach = approach,
                  adaptive = adaptive, kernel = kernel, p = p, theta = theta,
                  longlat = longlat, dMat = prep$dMat_sub)
}
