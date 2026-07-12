# Tests for the geographically weighted spatial lag (estimate_gwlag.R)

# Three collinear points at x = 0, 1, 2 (y = 0), values 10, 20, 30. With a boxcar
# kernel and a bandwidth larger than every distance, all in-range weights are 1,
# so the lag is a plain mean over the included neighbours - hand-computable.
.lag_fixture <- function() {
  data.frame(unit = c("a", "b", "c"),
             longitude = c(0, 1, 2), latitude = c(0, 0, 0),
             z = c(10, 20, 30), stringsAsFactors = FALSE)
}

test_that("neighbour mean excludes self by default (boxcar, flat weights)", {
  out <- estimate_gwlag_by_point(.lag_fixture(), unit = "unit", value_cols = "z",
                                 kernel = "boxcar", adaptive = FALSE, bw = 10,
                                 distance_metric = "Euclidean", include_self = FALSE)
  out <- out[match(c("a", "b", "c"), out$unit)]
  expect_s3_class(out, "data.table")
  expect_equal(out$z_LM, c(25, 20, 15))   # mean of the OTHER two values
})

test_that("include_self = TRUE averages over all units", {
  out <- estimate_gwlag_by_point(.lag_fixture(), unit = "unit", value_cols = "z",
                                 kernel = "boxcar", adaptive = FALSE, bw = 10,
                                 include_self = TRUE)
  out <- out[match(c("a", "b", "c"), out$unit)]
  expect_equal(out$z_LM, c(20, 20, 20))   # mean of all three
})

test_that("multiple value columns are lagged in one pass", {
  d <- .lag_fixture(); d$w <- c(1, 2, 3)
  out <- estimate_gwlag_by_point(d, unit = "unit", value_cols = c("z", "w"),
                                 kernel = "boxcar", adaptive = FALSE, bw = 10,
                                 include_self = FALSE)
  out <- out[match(c("a", "b", "c"), out$unit)]
  expect_true(all(c("z_LM", "w_LM") %in% names(out)))
  expect_equal(out$z_LM, c(25, 20, 15))
  expect_equal(out$w_LM, c(2.5, 2, 1.5))
  expect_equal(attr(out, "bandwidth"), 10)
})

test_that("by_polygon returns one lag row per polygon and finite values", {
  skip_if_not_installed("sf")
  g   <- make_sf_grid(n = 3L, id_col = "pid")     # 9 polygons p1..p9
  dat <- data.frame(pid = g$pid, z = seq_along(g$pid), stringsAsFactors = FALSE)
  out <- estimate_gwlag_by_polygon(dat, unit = "pid", polygons = g,
                                   value_cols = "z", kernel = "bisquare",
                                   adaptive = FALSE, bw = 5, include_self = FALSE)
  expect_s3_class(out, "data.table")
  expect_equal(sort(out$pid), sort(g$pid))
  expect_true(all(c("pid", "z_LM") %in% names(out)))
  expect_true(all(is.finite(out$z_LM)))
})

test_that("by_polygon errors on a non-polygon input", {
  expect_error(
    estimate_gwlag_by_polygon(data.frame(pid = "p1", z = 1), unit = "pid",
                              polygons = "nope", value_cols = "z"))
})
