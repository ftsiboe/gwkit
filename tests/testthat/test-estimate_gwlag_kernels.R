# estimate_gwlag_kernels(): the multi-kernel path must reproduce looping
# estimate_gwlag() kernel-by-kernel (per_kernel), reusing one distance build;
# "shared" reuses a single bandwidth; a named bw is used verbatim.

.lagk_points <- function(n = 5L) {
  g <- expand.grid(longitude = seq_len(n), latitude = seq_len(n))
  data.frame(unit = paste0("u", seq_len(nrow(g))),
             longitude = g$longitude, latitude = g$latitude,
             z = as.numeric(seq_len(nrow(g))), stringsAsFactors = FALSE)
}

test_that("kernel column + per-kernel bandwidth attribute are returned", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .lagk_points()
  ks <- c("gaussian", "bisquare", "boxcar")
  out <- estimate_gwlag_kernels(d, unit = "unit", value_cols = "z",
    kernel = ks, adaptive = FALSE, bw = 3, distance_metric = "Euclidean",
    include_self = FALSE)
  expect_s3_class(out, "data.table")
  expect_true(all(c("kernel", "z_LM") %in% names(out)))
  expect_setequal(unique(out$kernel), ks)
  expect_length(attr(out, "bandwidth"), length(ks))
  expect_setequal(names(attr(out, "bandwidth")), ks)
})

test_that("per_kernel matches looping estimate_gwlag", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .lagk_points()
  ks <- c("gaussian", "bisquare")
  batch <- estimate_gwlag_kernels(d, unit = "unit", value_cols = "z",
    kernel = ks, adaptive = FALSE, bw = NULL, bandwidth = "per_kernel",
    distance_metric = "Euclidean", include_self = FALSE)
  for (k in ks) {
    single <- data.table::as.data.table(estimate_gwlag(d, unit = "unit", value_cols = "z",
      kernel = k, adaptive = FALSE, bw = NULL, distance_metric = "Euclidean",
      include_self = FALSE))
    b <- data.table::copy(batch[kernel == k])
    data.table::setorder(b, unit); data.table::setorder(single, unit)
    expect_equal(b$z_LM, single$z_LM, tolerance = 1e-9)
    expect_equal(unname(attr(batch, "bandwidth")[k]),
                 unname(attr(single, "bandwidth")), tolerance = 1e-9)
  }
})

test_that("shared reuses one bandwidth across kernels", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  out <- estimate_gwlag_kernels(.lagk_points(), unit = "unit", value_cols = "z",
    kernel = c("gaussian", "bisquare", "tricube"), adaptive = FALSE, bw = NULL,
    bandwidth = "shared", distance_metric = "Euclidean", include_self = FALSE)
  expect_length(unique(unname(attr(out, "bandwidth"))), 1L)
})

test_that("flat weights recover the self-excluded neighbour mean", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  out <- estimate_gwlag_kernels(
    data.frame(unit = c("a","b","c"), longitude = c(0,1,2), latitude = c(0,0,0), z = c(10,20,30)),
    unit = "unit", value_cols = "z", kernel = "boxcar", adaptive = FALSE, bw = 1e6,
    include_self = FALSE)
  b <- out[kernel == "boxcar"][match(c("a","b","c"), unit)]
  expect_equal(b$z_LM, c(25, 20, 15), tolerance = 1e-9)
})

test_that("a named bw is used verbatim and round-trips", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .lagk_points(); ks <- c("gaussian", "bisquare")
  first  <- estimate_gwlag_kernels(d, unit = "unit", value_cols = "z", kernel = ks,
    adaptive = FALSE, bw = NULL, bandwidth = "per_kernel", include_self = FALSE)
  bw_vec <- attr(first, "bandwidth")
  second <- estimate_gwlag_kernels(d, unit = "unit", value_cols = "z", kernel = ks,
    adaptive = FALSE, bw = bw_vec, include_self = FALSE)
  data.table::setorder(first, kernel, unit); data.table::setorder(second, kernel, unit)
  expect_equal(second$z_LM, first$z_LM, tolerance = 1e-12)
  expect_equal(attr(second, "bandwidth"), bw_vec)
})
