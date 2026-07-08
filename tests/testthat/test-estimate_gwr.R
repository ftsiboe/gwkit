# Uses the real spatial layout of GWmodel's LondonHP points, extended into a
# short per-unit panel so there is a covariate to regress the response on.

load_londonhp <- function() {
  e <- new.env()
  utils::data("LondonHP", package = "GWmodel", envir = e)
  e$londonhp
}

test_that("estimate_gwr_by_point returns the expected schema and finite fits", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  sp_pts <- load_londonhp()
  co     <- sp::coordinates(sp_pts)

  set.seed(1L)
  n   <- min(60L, nrow(co))
  sel <- sample.int(nrow(co), n)
  b   <- seq_len(n) / n                       # a distinct "true" slope per unit
  tt  <- 0:5                                  # covariate values within each unit

  panel <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      unit  = paste0("u", i),
      lon   = co[sel[i], 1],
      lat   = co[sel[i], 2],
      trend = tt,
      y     = 10 + b[i] * tt + rnorm(length(tt), sd = 0.05),
      stringsAsFactors = FALSE
    )
  }))

  res <- estimate_gwr_by_point(
    panel, unit = "unit", response = "y", covariate = "trend",
    coords = c("lon", "lat"), distance_metric = "Euclidean",
    kernel = "bisquare", adaptive = FALSE, variance = TRUE
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), n)
  expect_true(all(c("term", "unit_id", "unit_level",
                    "mean_estimate", "mean_standard_error", "mean_t_value",
                    "mean_p_value", "var_estimate", "var_p_value",
                    "model_estimator") %in% names(res)))
  expect_identical(sort(res$unit_id), sort(paste0("u", seq_len(n))))
  expect_true(all(res$model_estimator == "gwr"))
  expect_true(all(res$unit_level == "unit"))

  # local slopes should be recovered for the bulk of units (allow a few NA at
  # degenerate/edge locations)
  expect_gt(mean(is.finite(res$mean_estimate)), 0.5)
  expect_gt(mean(is.finite(res$mean_p_value)),  0.5)
})

test_that("variance = FALSE omits the variance columns", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  sp_pts <- load_londonhp()
  co     <- sp::coordinates(sp_pts)

  set.seed(2L)
  n   <- min(40L, nrow(co))
  sel <- sample.int(nrow(co), n)
  tt  <- 0:4
  panel <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(unit = paste0("u", i), lon = co[sel[i], 1], lat = co[sel[i], 2],
               trend = tt, y = 5 + 0.2 * tt + rnorm(length(tt), sd = 0.05))
  }))

  res <- estimate_gwr_by_point(
    panel, unit = "unit", response = "y", covariate = "trend",
    coords = c("lon", "lat"), distance_metric = "Euclidean",
    adaptive = FALSE, variance = FALSE
  )
  expect_true("mean_estimate" %in% names(res))
  expect_false("var_estimate" %in% names(res))
})
