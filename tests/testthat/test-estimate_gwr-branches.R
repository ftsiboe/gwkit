# Exercises the estimate_gwr branches the LondonHP test does not reach:
# every kernel, the adaptive bandwidth path, the supplied-bw path, the unknown
# kernel error, the too-few-observations guard, and the polygon variant.

test_that("every kernel runs and yields a slope per unit (supplied fixed bw)", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  panel <- make_point_panel()
  for (k in c("gaussian", "exponential", "bisquare", "boxcar", "tricube")) {
    res <- estimate_gwr_by_point(
      panel, unit = "unit", response = "y", covariate = "trend",
      coords = c("lon", "lat"), kernel = k, adaptive = FALSE,
      distance_metric = "Euclidean", bw = 3, variance = FALSE
    )
    expect_s3_class(res, "data.table")
    expect_equal(nrow(res), 25L)
    expect_gt(mean(is.finite(res$mean_estimate)), 0.5)
  }
})

test_that("adaptive bandwidth path runs and variance = TRUE adds var_ columns", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  panel <- make_point_panel()
  res <- estimate_gwr_by_point(
    panel, unit = "unit", response = "y", covariate = "trend",
    coords = c("lon", "lat"), kernel = "bisquare", adaptive = TRUE,
    distance_metric = "Euclidean", bw = 30, variance = TRUE
  )
  expect_equal(nrow(res), 25L)
  expect_true(all(c("var_estimate", "var_standard_error", "var_p_value") %in% names(res)))
})

test_that("an unknown kernel is rejected", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  panel <- make_point_panel()
  expect_error(
    estimate_gwr_by_point(
      panel, unit = "unit", response = "y", covariate = "trend",
      coords = c("lon", "lat"), kernel = "nope", adaptive = FALSE,
      distance_metric = "Euclidean", bw = 3
    ),
    "Unknown kernel"
  )
})

test_that("too few finite observations returns NULL", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  small <- data.frame(unit = "u1", lon = 1, lat = 1, trend = 0:3, y = c(1, 2, 3, 4))
  expect_null(
    estimate_gwr_by_point(
      small, unit = "unit", response = "y", covariate = "trend",
      coords = c("lon", "lat"), distance_metric = "Euclidean", bw = 1
    )
  )
})

test_that("estimate_gwr_by_polygon fits at polygon centroids", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  skip_if_not_installed("sf")

  pg  <- make_sf_grid(5L, id_col = "pid")
  ids <- pg$pid
  tt  <- 0:4
  set.seed(2L)
  dat <- do.call(rbind, lapply(seq_along(ids), function(i) {
    data.frame(
      pid = ids[i], trend = tt,
      y   = 5 + 0.1 * i * tt + stats::rnorm(length(tt), sd = 0.05),
      stringsAsFactors = FALSE
    )
  }))

  res <- estimate_gwr_by_polygon(
    dat, unit = "pid", polygons = pg, response = "y", covariate = "trend",
    poly_id = "pid", distance_metric = "Euclidean", kernel = "bisquare",
    adaptive = FALSE, bw = 3, variance = TRUE
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), length(ids))
  expect_true(all(c("mean_estimate", "var_estimate", "model_estimator") %in% names(res)))
  expect_identical(sort(res$unit_id), sort(ids))
})

test_that("estimate_gwr_by_polygon rejects a non-spatial polygons argument", {
  expect_error(
    estimate_gwr_by_polygon(
      data.frame(pid = "p1", trend = 0:4, y = 1:5),
      unit = "pid", polygons = data.frame(pid = "p1"), response = "y"
    ),
    "SpatVector or an sf"
  )
})
