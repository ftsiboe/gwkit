# Exercises the estimate_gwss_by_polygon body (join, points-on-surface, bandwidth
# CV, gwss) on a small synthetic sf lattice, plus its "too few observed" guard.
#
# "Great Circle" is used so the internal reprojection is a no-op to EPSG:4326.
# With a non-longlat metric the function projects to `target_crs` (5070, a US
# CRS), which would distort this synthetic lat/lon lattice.

test_that("estimate_gwss_by_polygon returns local summaries at all polygons", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")
  skip_if_not_installed("dplyr")

  pg <- make_sf_grid(5L, id_col = "pid")
  set.seed(3L)
  dat <- data.frame(pid = pg$pid, v = stats::rnorm(nrow(pg), mean = 10),
                    stringsAsFactors = FALSE)

  res <- estimate_gwss_by_polygon(
    data            = dat,
    shape_file      = pg,
    fip_col         = "pid",
    variable_list   = "v",
    distance_metric = "Great Circle",
    kernel          = "gaussian",
    adaptive        = TRUE
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), nrow(pg))
  expect_true("pid" %in% names(res))
  # GWmodel names local means <var>_LM
  expect_true(any(grepl("_LM$", names(res))))
})

test_that("estimate_gwss_by_polygon returns NULL when too few polygons are observed", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dplyr")

  pg  <- make_sf_grid(3L, id_col = "pid")
  dat <- data.frame(pid = pg$pid[1:3], v = c(1, 2, 3), stringsAsFactors = FALSE)

  expect_message(
    res <- estimate_gwss_by_polygon(dat, pg, "pid", "v",
                                    distance_metric = "Great Circle"),
    "Too few"
  )
  expect_null(res)
})
