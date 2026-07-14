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
