# estimate_gwss() input validation - hits the stop() branches without touching
# GWmodel, so these are fast and dependency-light.

test_that("estimate_gwss validates shared inputs", {
  df <- data.frame(lon = as.numeric(1:10), lat = as.numeric(1:10),
                   value = stats::rnorm(10))

  expect_error(estimate_gwss(), "must be supplied")                        # missing data
  expect_error(estimate_gwss(df, "value", kernel = "bad"),   "kernel")
  expect_error(estimate_gwss(df, "value", approach = "bad"), "approach")

  # point mode: variable column missing from data
  expect_error(estimate_gwss(df, "missing", coords = c("lon", "lat")),
               "must contain")
  # a non-spatial geometry is rejected
  expect_error(estimate_gwss(df, "value", coords = c("lon", "lat"),
                             geometry = data.frame(x = 1)),
               "SpatVector or an sf")
})

test_that("estimate_gwss validates polygon-mode inputs", {
  skip_if_not_installed("sf")

  pg  <- make_sf_grid(3L, id_col = "pid")
  dat <- data.frame(pid = pg$pid, v = stats::rnorm(nrow(pg)), stringsAsFactors = FALSE)

  # id_col is required when a geometry is supplied
  expect_error(estimate_gwss(dat, "v", geometry = pg), "id_col")
  # data must contain the variable column
  expect_error(estimate_gwss(dat, "missing", geometry = pg, id_col = "pid"),
               "must contain")
})
