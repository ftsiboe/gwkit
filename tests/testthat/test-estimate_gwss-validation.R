# The three GWSS entry points share a long validation preamble. These tests hit
# each stop() branch without touching GWmodel, so they are fast and dependency-free.

test_that("estimate_gwss_by_point validates its inputs", {
  df <- data.frame(lon = as.numeric(1:10), lat = as.numeric(1:10),
                   value = stats::rnorm(10))

  expect_error(estimate_gwss_by_point(), "locat_observed")
  expect_error(estimate_gwss_by_point(df), "locat_predict")
  expect_error(estimate_gwss_by_point(df, df, NULL, "lat", "value"),
               "character column names")
  expect_error(estimate_gwss_by_point(df, df, "lon", "lat", "value", kernel = "bad"),
               "kernel")
  expect_error(estimate_gwss_by_point(df, df, "lon", "lat", "value", approach = "bad"),
               "approach")
  expect_error(estimate_gwss_by_point(df, df, "lon", "lat", "missing"),
               "must contain")
  expect_error(
    estimate_gwss_by_point(df, df, "lon", "lat", "value", group_variable = "nope"),
    "group_variable"
  )
  expect_error(
    estimate_gwss_by_point(df, df, "lon", "lat", "value", identifiers = "nope"),
    "identifiers"
  )
})

test_that("estimate_gwss_by_polygon validates its inputs", {
  skip_if_not_installed("sf")

  pg  <- make_sf_grid(3L, id_col = "pid")
  dat <- data.frame(pid = pg$pid, v = stats::rnorm(nrow(pg)),
                    stringsAsFactors = FALSE)

  expect_error(estimate_gwss_by_polygon(), "`data` must be supplied")
  expect_error(estimate_gwss_by_polygon(dat), "`shape_file` must be supplied")
  expect_error(estimate_gwss_by_polygon(dat, pg, NULL, "v"), "missing one of")
  expect_error(estimate_gwss_by_polygon(dat, pg, "pid", "v", kernel = "bad"), "kernel")
  expect_error(estimate_gwss_by_polygon(dat, pg, "pid", "v", approach = "bad"), "approach")
  expect_error(estimate_gwss_by_polygon(dat, pg, "pid", "missing"), "must contain columns")
  expect_error(estimate_gwss_by_polygon(dat, data.frame(pid = "p1"), "pid", "v"),
               "must be an `sf` object")
})

test_that("estimate_gwss_by_county validates its inputs", {
  dat <- data.frame(fips = c("01001", "01003"), v = c(1, 2), stringsAsFactors = FALSE)

  expect_error(estimate_gwss_by_county(), "`data` must be supplied")
  expect_error(estimate_gwss_by_county(dat, NULL, "v"), "missing one of")
  expect_error(estimate_gwss_by_county(dat, "fips", "v", kernel = "bad"), "kernel")
  expect_error(estimate_gwss_by_county(dat, "fips", "v", approach = "bad"), "approach")
  expect_error(estimate_gwss_by_county(dat, "fips", "missing"), "must contain columns")
})
