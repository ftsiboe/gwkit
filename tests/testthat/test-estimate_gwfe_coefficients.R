# Tests for GW fixed-effects local regression (estimate_gwfe_coefficients.R)

# A panel: 6 units over 5 years, y_it = alpha_i + 2 x_it + tiny noise. The within
# (fixed-effects) slope on x is ~2. With a boxcar kernel and a huge bandwidth all
# weights are 1, so every focal unit recovers the SAME global within slope.
.fe_panel <- function(seed = 1L) {
  set.seed(seed)
  units <- paste0("u", 1:6)
  do.call(rbind, lapply(seq_along(units), function(i) {
    x <- rnorm(5)
    data.frame(unit = units[i], longitude = i, latitude = 0, year = 1:5,
               x = x, y = i * 10 + 2 * x + rnorm(5, sd = 0.01),
               stringsAsFactors = FALSE)
  }))
}

# reference within slope via manual demeaning
.within_slope <- function(d) {
  d$x_dm <- d$x - ave(d$x, d$unit)
  d$y_dm <- d$y - ave(d$y, d$unit)
  unname(coef(lm(y_dm ~ x_dm, data = d))["x_dm"])
}

test_that("flat weights recover the national within (FE) slope at every unit", {
  d  <- .fe_panel()
  bw_ref <- .within_slope(d)

  out <- estimate_gwfe_coefficients_by_point(
    d, unit = "unit", formula = y ~ x, panel = "unit", time = "year",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, fe_df_correction = FALSE)

  expect_s3_class(out, "data.table")
  xs <- out[term == "x", est]
  expect_equal(length(xs), length(unique(d$unit)))
  expect_equal(xs, rep(bw_ref, length(xs)), tolerance = 1e-6)
  expect_true(all(out$model_estimator == "gwfe"))
  expect_equal(unname(attr(out, "NFE")), length(unique(d$unit)))
})

test_that("local within-R^2 and fit columns are present and finite", {
  d <- .fe_panel()
  out <- estimate_gwfe_coefficients_by_point(
    d, unit = "unit", formula = y ~ x, panel = "unit", time = "year",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, terms = "x")
  expect_true(all(c("r_squared","adj_r_squared","n_obs","n_units") %in% names(out)))
  expect_true(all(is.finite(out$r_squared)))
  expect_true(all(out$r_squared >= 0 & out$r_squared <= 1))
  expect_true(all(out$n_units >= 1))
})

test_that("FE df correction inflates the SE relative to the naive lm df", {
  d <- .fe_panel()
  raw <- estimate_gwfe_coefficients_by_point(
    d, unit = "unit", formula = y ~ x, panel = "unit", time = "year",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, terms = "x",
    fe_df_correction = FALSE)
  adj <- estimate_gwfe_coefficients_by_point(
    d, unit = "unit", formula = y ~ x, panel = "unit", time = "year",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, terms = "x",
    fe_df_correction = TRUE)
  expect_true(all(adj[term == "x", se] >= raw[term == "x", se] - 1e-12))
})

test_that("multi-covariate formula returns a slope per covariate", {
  set.seed(3); d <- .fe_panel()
  d$z <- rnorm(nrow(d))
  out <- estimate_gwfe_coefficients_by_point(
    d, unit = "unit", formula = y ~ x + z, panel = "unit", time = "year",
    kernel = "boxcar", adaptive = FALSE, bw = 1e6)
  expect_setequal(unique(out$term), c("(Intercept)", "x", "z"))
})

test_that("by_polygon runs over an sf grid", {
  skip_if_not_installed("sf")
  g <- make_sf_grid(n = 3L, id_col = "pid")     # 9 polygons
  set.seed(4)
  d <- do.call(rbind, lapply(g$pid, function(id) {
    x <- rnorm(5)
    data.frame(pid = id, year = 1:5, x = x, y = rnorm(1)*10 + 1.5 * x + rnorm(5, sd = 0.05),
               stringsAsFactors = FALSE)
  }))
  out <- estimate_gwfe_coefficients_by_polygon(
    d, unit = "pid", polygons = g, formula = y ~ x, panel = "pid", time = "year",
    kernel = "bisquare", adaptive = TRUE, bw = 50, terms = "x")
  expect_s3_class(out, "data.table")
  expect_true(all(c("term","unit_id","est","r_squared") %in% names(out)))
  expect_true("x" %in% out$term)
})
