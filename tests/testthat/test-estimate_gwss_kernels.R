# estimate_gwss_kernels(): the cached multi-kernel path must reproduce looping
# estimate_gwss() kernel-by-kernel (bandwidth = "per_kernel"), since it reuses the
# same seeded subsample and the same distance matrices. "shared" reuses a single
# bandwidth across the kernels.

.londonhp_points <- function() {
  e <- new.env()
  utils::data("LondonHP", package = "GWmodel", envir = e)
  obs <- sf::st_as_sf(e$londonhp)
  suppressWarnings(sf::st_crs(obs) <- 27700)   # LondonHP is on the OSGB grid
  obs <- sf::st_transform(obs, 4326)
  co  <- sf::st_coordinates(obs)
  data.frame(lon = co[, 1], lat = co[, 2], value = as.numeric(obs$PURCHASE))
}

test_that("estimate_gwss_kernels (per_kernel) matches looping estimate_gwss", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  df      <- .londonhp_points()
  kernels <- c("gaussian", "bisquare")
  args <- list(data = df, variable_list = "value", coords = c("lon", "lat"),
               distance_metric = "Great Circle", adaptive = TRUE,
               target_crs = 27700, feasible_radius = 500)

  batch <- do.call(estimate_gwss_kernels,
                   c(args, list(kernel = kernels, bandwidth = "per_kernel")))
  expect_s3_class(batch, "data.table")
  skip_if(nrow(batch) == 0, "GWSS produced no rows in this environment.")
  expect_true("kernel" %in% names(batch))
  expect_setequal(unique(batch$kernel), kernels)

  for (k in kernels) {
    single <- data.table::as.data.table(do.call(estimate_gwss, c(args, list(kernel = k))))
    b <- data.table::copy(batch[kernel == k])
    data.table::setorder(b, longitude, latitude)
    data.table::setorder(single, longitude, latitude)
    expect_equal(nrow(b), nrow(single))
    # Local mean surface must be identical to the per-kernel single call.
    expect_equal(b$value_LM, single$value_LM, tolerance = 1e-9)
  }
})

test_that("estimate_gwss_kernels (shared) reuses one bandwidth across kernels", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  df  <- .londonhp_points()
  res <- estimate_gwss_kernels(
    data = df, variable_list = "value", coords = c("lon", "lat"),
    distance_metric = "Great Circle", kernel = c("gaussian", "bisquare"),
    adaptive = TRUE, target_crs = 27700, feasible_radius = 500,
    bandwidth = "shared")

  expect_s3_class(res, "data.table")
  skip_if(nrow(res) == 0, "GWSS produced no rows in this environment.")
  expect_true("kernel" %in% names(res))
  # A single shared bandwidth is attached as an attribute.
  expect_false(is.null(attr(res, "bandwidth")))
  expect_length(attr(res, "bandwidth"), 1L)
})
