# estimate_gwr_kernels(): the cached multi-kernel path must reproduce looping
# estimate_gwr() kernel-by-kernel (it reuses the same distance build and, under
# bandwidth = "per_kernel", the same per-kernel bw.gwr selection). "shared"
# reuses a single bandwidth across the kernels. A named `bw` is fed back verbatim.

# A panel: n units on a line (lon = i, lat = 0), each observed 5 times with a
# unit fixed effect and a common within slope, so FE GWR is well-posed.
.gwrk_fe_panel <- function(n = 8L, seed = 1L) {
  set.seed(seed)
  units <- paste0("u", seq_len(n))
  do.call(rbind, lapply(seq_along(units), function(i) {
    x <- rnorm(5); z <- rnorm(5)
    data.frame(unit = units[i], lon = i, lat = 0, year = 1:5,
               x = x, z = z, y = i * 10 + 2 * x - 1.5 * z + rnorm(5, sd = 0.01),
               stringsAsFactors = FALSE)
  }))
}

test_that("kernel column + per-kernel bandwidth attribute are returned", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwrk_fe_panel()
  ks <- c("gaussian", "bisquare", "boxcar")
  out <- estimate_gwr_kernels(d, unit = "unit", formula = y ~ x + z,
    panel = "unit", time = "year", coords = c("lon", "lat"),
    kernel = ks, adaptive = TRUE, bandwidth = "per_kernel", fit_stats = TRUE)

  expect_s3_class(out, "data.table")
  expect_true("kernel" %in% names(out))
  expect_setequal(unique(out$kernel), ks)
  expect_true(all(out$estimand == "mean"))
  expect_true(all(out$model_estimator == "gwfe"))
  bwv <- attr(out, "bandwidth")
  expect_length(bwv, length(ks))
  expect_setequal(names(bwv), ks)
})

test_that("per_kernel matches looping estimate_gwr kernel-by-kernel", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .gwrk_fe_panel()
  ks <- c("gaussian", "bisquare")
  args <- list(data = d, unit = "unit", formula = y ~ x + z,
               panel = "unit", time = "year", coords = c("lon", "lat"),
               adaptive = TRUE, fit_stats = TRUE)

  batch <- do.call(estimate_gwr_kernels,
                   c(args, list(kernel = ks, bandwidth = "per_kernel")))

  for (k in ks) {
    single <- data.table::as.data.table(do.call(estimate_gwr, c(args, list(kernel = k))))
    b <- data.table::copy(batch[kernel == k])
    data.table::setorder(b, unit_id, term)
    data.table::setorder(single, unit_id, term)
    expect_equal(nrow(b), nrow(single))
    expect_equal(b$estimate,  single$estimate,  tolerance = 1e-9)
    expect_equal(b$r_squared, single$r_squared, tolerance = 1e-9)
    # the selected bandwidth must match the single call's, too
    expect_equal(unname(attr(batch, "bandwidth")[k]),
                 unname(attr(single, "bandwidth")), tolerance = 1e-9)
  }
})

test_that("shared reuses one bandwidth across kernels", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .gwrk_fe_panel()
  out <- estimate_gwr_kernels(d, unit = "unit", formula = y ~ x + z,
    panel = "unit", time = "year", coords = c("lon", "lat"),
    kernel = c("gaussian", "bisquare", "tricube"), adaptive = TRUE,
    bandwidth = "shared")
  bwv <- attr(out, "bandwidth")
  expect_length(unique(unname(bwv)), 1L)   # one value, reused
})

test_that("flat weights recover the global within slope at every unit & kernel", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwrk_fe_panel()
  d$x_dm <- d$x - ave(d$x, d$unit); d$z_dm <- d$z - ave(d$z, d$unit)
  d$y_dm <- d$y - ave(d$y, d$unit)
  within <- coef(lm(y_dm ~ x_dm + z_dm, data = d))

  out <- estimate_gwr_kernels(d, unit = "unit", formula = y ~ x + z,
    panel = "unit", time = "year", coords = c("lon", "lat"),
    kernel = c("boxcar", "gaussian"), adaptive = FALSE, bw = 1e6, terms = c("x", "z"))

  for (k in unique(out$kernel)) {
    xs <- out[kernel == k & term == "x", estimate]
    zs <- out[kernel == k & term == "z", estimate]
    expect_equal(xs, rep(unname(within["x_dm"]), length(xs)), tolerance = 1e-6)
    expect_equal(zs, rep(unname(within["z_dm"]), length(zs)), tolerance = 1e-6)
  }
})

test_that("a named bw is used verbatim (no selection) and round-trips", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d  <- .gwrk_fe_panel()
  ks <- c("gaussian", "bisquare")
  args <- list(data = d, unit = "unit", formula = y ~ x + z,
               panel = "unit", time = "year", coords = c("lon", "lat"),
               adaptive = TRUE, fit_stats = TRUE, kernel = ks)

  first  <- do.call(estimate_gwr_kernels, c(args, list(bandwidth = "per_kernel")))
  bw_vec <- attr(first, "bandwidth")                     # named per-kernel bws
  second <- do.call(estimate_gwr_kernels, c(args, list(bw = bw_vec)))

  data.table::setorder(first,  kernel, unit_id, term)
  data.table::setorder(second, kernel, unit_id, term)
  expect_equal(second$estimate, first$estimate, tolerance = 1e-12)
  expect_equal(attr(second, "bandwidth"), bw_vec)
})
