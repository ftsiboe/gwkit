# Tests for multi-covariate local GWR coefficients (estimate_gwr_coefficients.R)

# With a boxcar kernel and a bandwidth larger than every distance, all weights
# are 1, so every focal unit fits the SAME ordinary lm over all observations.
# The local coefficients must therefore equal the global OLS fit at every unit.
.coef_points <- function(n = 6L, seed = 1L) {
  set.seed(seed)
  data.frame(unit = paste0("u", seq_len(n)),
             longitude = seq_len(n), latitude = 0,
             x1 = rnorm(n), x2 = rnorm(n),
             y  = rnorm(n), stringsAsFactors = FALSE)
}

test_that("flat weights reproduce the global OLS coefficients at every unit", {
  d   <- .coef_points()
  ols <- coef(lm(y ~ x1 + x2, data = d))

  out <- estimate_gwr_coefficients_by_point(
    d, unit = "unit", formula = y ~ x1 + x2,
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, distance_metric = "Euclidean")

  expect_s3_class(out, "data.table")
  expect_setequal(unique(out$term), c("(Intercept)", "x1", "x2"))
  expect_equal(sort(unique(out$unit_id)), sort(d$unit))

  for (tm in names(ols)) {
    vals <- out[term == tm, est]
    expect_equal(vals, rep(unname(ols[tm]), length(vals)), tolerance = 1e-8)
  }
  expect_true(all(out$model_estimator == "gwr_wls"))
})

test_that("single-covariate slope matches lm, and terms filter works", {
  d   <- .coef_points()
  ols <- coef(lm(y ~ x1, data = d))["x1"]

  out <- estimate_gwr_coefficients_by_point(
    d, unit = "unit", formula = y ~ x1, terms = "x1",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6)

  expect_equal(unique(out$term), "x1")
  expect_equal(out$est, rep(unname(ols), nrow(out)), tolerance = 1e-8)
})

test_that("standard errors match stats::lm exactly under flat weights", {
  d  <- .coef_points()
  sm <- coef(summary(lm(y ~ x1 + x2, data = d)))
  out <- estimate_gwr_coefficients_by_point(
    d, unit = "unit", formula = y ~ x1 + x2,
    kernel = "boxcar", adaptive = FALSE, bw = 1e6)
  one <- out[unit_id == "u1"][match(rownames(sm), term)]
  expect_equal(one$se, unname(sm[, "Std. Error"]), tolerance = 1e-8)
})

test_that("by_polygon returns a term x unit table over an sf grid", {
  skip_if_not_installed("sf")
  g <- make_sf_grid(n = 3L, id_col = "pid")
  set.seed(2)
  dat <- data.frame(pid = g$pid, x = rnorm(length(g$pid)),
                    y = rnorm(length(g$pid)), stringsAsFactors = FALSE)
  out <- estimate_gwr_coefficients_by_polygon(
    dat, unit = "pid", polygons = g, formula = y ~ x,
    kernel = "bisquare", adaptive = TRUE, bw = 50)
  expect_s3_class(out, "data.table")
  expect_true(all(c("term", "unit_id", "est", "se", "tv", "pv") %in% names(out)))
  expect_true("x" %in% out$term)
})
