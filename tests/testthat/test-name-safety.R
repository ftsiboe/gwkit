# Robustness regressions:
#  (1) internal coordinate/id columns are dot-prefixed, so user model columns
#      literally named X / Y must not collide with them (polygon mode).
#  (2) estimate_gwss() must not disturb the caller's RNG stream.

test_that("model columns named X and Y do not collide with internal coords (polygon)", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sf"); skip_if_not_installed("sp")

  pg  <- make_sf_grid(5L, id_col = "pid")
  set.seed(7L)
  dat <- data.frame(pid = pg$pid,
                    Y = stats::rnorm(nrow(pg)),     # response literally named "Y"
                    X = stats::rnorm(nrow(pg)),     # covariate literally named "X"
                    stringsAsFactors = FALSE)

  out <- estimate_gwr(dat, unit = "pid", formula = Y ~ X,
                      geometry = pg, poly_id = "pid",
                      distance_metric = "Great Circle",
                      kernel = "boxcar", adaptive = FALSE, bw = 1e6)

  expect_s3_class(out, "data.table")
  expect_setequal(unique(out$term), c("(Intercept)", "X"))
  expect_equal(nrow(out[estimand == "mean" & term == "X"]), nrow(pg))
  expect_true(all(is.finite(out$estimate)))
})

test_that("estimate_gwss does not disturb the caller's RNG stream", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sf"); skip_if_not_installed("sp")

  pg  <- make_sf_grid(5L, id_col = "pid")
  set.seed(3L)
  dat <- data.frame(pid = pg$pid,
                    a = stats::rnorm(nrow(pg), mean = 10),
                    b = stats::rnorm(nrow(pg), mean = 10),
                    stringsAsFactors = FALSE)

  set.seed(123L)
  before <- .Random.seed
  invisible(estimate_gwss(dat, variable_list = c("a", "b"), geometry = pg,
                          id_col = "pid", distance_metric = "Great Circle",
                          kernel = "gaussian", adaptive = TRUE))
  expect_identical(.Random.seed, before)
})
