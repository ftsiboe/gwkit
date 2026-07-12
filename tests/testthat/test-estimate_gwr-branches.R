# estimate_gwr(): kernels, adaptive, variance mode, errors, polygon (+FE) variants.

test_that("every kernel runs and yields a coefficient per unit (supplied bw)", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  panel <- make_point_panel()
  for (k in c("gaussian", "exponential", "bisquare", "boxcar", "tricube")) {
    res <- estimate_gwr(panel, unit = "unit", formula = y ~ trend,
      coords = c("lon", "lat"), kernel = k, adaptive = FALSE, bw = 3, terms = "trend")
    expect_s3_class(res, "data.table")
    expect_equal(nrow(res), 25L)
    expect_gt(mean(is.finite(res$estimate)), 0.5)
  }
})

test_that("adaptive bandwidth path runs", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  panel <- make_point_panel()
  res <- estimate_gwr(panel, unit = "unit", formula = y ~ trend,
    coords = c("lon", "lat"), kernel = "bisquare", adaptive = TRUE, bw = 30, terms = "trend")
  expect_equal(nrow(res), 25L)
})

test_that("variance = TRUE adds estimand == 'variance' rows", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  panel <- make_point_panel()
  res <- estimate_gwr(panel, unit = "unit", formula = y ~ trend,
    coords = c("lon", "lat"), kernel = "bisquare", adaptive = FALSE, bw = 3,
    variance = TRUE, terms = "trend")
  expect_setequal(unique(res$estimand), c("mean", "variance"))
  expect_equal(nrow(res[estimand == "mean"]), 25L)
  expect_equal(nrow(res[estimand == "variance"]), 25L)
})

test_that("an unknown kernel is rejected", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  panel <- make_point_panel()
  expect_error(
    estimate_gwr(panel, unit = "unit", formula = y ~ trend,
      coords = c("lon", "lat"), kernel = "nope", adaptive = FALSE, bw = 3),
    "Unknown kernel")
})

test_that("too few finite observations returns NULL", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  small <- data.frame(unit = "u1", lon = 1, lat = 1, trend = 0:1, y = c(1, 2))
  expect_null(
    estimate_gwr(small, unit = "unit", formula = y ~ trend,
      coords = c("lon", "lat"), bw = 1, distance_metric = "Euclidean"))
})

test_that("by_polygon runs plain and FE over an sf grid", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp"); skip_if_not_installed("sf")
  g   <- make_sf_grid(5L, id_col = "pid")
  ids <- g$pid; tt <- 0:4; set.seed(2)
  dat <- do.call(rbind, lapply(seq_along(ids), function(i) {
    data.frame(pid = ids[i], year = tt, x = rnorm(length(tt)),
               y = 5 + 0.1 * i * tt + rnorm(length(tt), sd = 0.05),
               stringsAsFactors = FALSE)
  }))

  res <- estimate_gwr(dat, unit = "pid", geometry =g, formula = y ~ x,
    kernel = "bisquare", adaptive = FALSE, bw = 3, terms = "x")
  expect_s3_class(res, "data.table")
  expect_identical(sort(unique(res$unit_id)), sort(ids))
  expect_true(all(res$model_estimator == "gwr"))

  resfe <- estimate_gwr(dat, unit = "pid", geometry =g, formula = y ~ x,
    panel = "pid", time = "year", kernel = "bisquare", adaptive = TRUE, bw = 20,
    fit_stats = TRUE, terms = "x")
  expect_true(all(resfe$model_estimator == "gwfe"))
  expect_true("r_squared" %in% names(resfe))
})

test_that("by_polygon rejects a non-spatial polygons argument", {
  expect_error(
    estimate_gwr(
      data.frame(pid = "p1", trend = 0:4, y = 1:5),
      unit = "pid", geometry =data.frame(pid = "p1"), formula = y ~ trend),
    "SpatVector or an sf")
})
