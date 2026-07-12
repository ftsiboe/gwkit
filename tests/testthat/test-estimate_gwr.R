# Unified estimate_gwr(): plain GWR, fixed-effects GWR, fit stats, opt-in SEs.

# A panel: n units (lon = i, lat = 0), each observed 5 times with y = a_i + 2x.
.gwr_fe_panel <- function(n = 6L, seed = 1L) {
  set.seed(seed)
  units <- paste0("u", seq_len(n))
  do.call(rbind, lapply(seq_along(units), function(i) {
    x <- rnorm(5)
    data.frame(unit = units[i], lon = i, lat = 0, year = 1:5,
               x = x, y = i * 10 + 2 * x + rnorm(5, sd = 0.01),
               stringsAsFactors = FALSE)
  }))
}

test_that("plain GWR: flat weights recover global OLS at every unit (estimand = mean)", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwr_fe_panel(); d$z <- rnorm(nrow(d))
  ols <- coef(lm(y ~ x + z, data = d))

  out <- estimate_gwr(d, unit = "unit", formula = y ~ x + z,
    coords = c("lon", "lat"), kernel = "boxcar", adaptive = FALSE, bw = 1e6)

  expect_s3_class(out, "data.table")
  expect_true(all(out$estimand == "mean"))
  expect_true(all(out$model_estimator == "gwr"))
  expect_setequal(unique(out$term), c("(Intercept)", "x", "z"))
  for (tm in names(ols)) {
    v <- out[term == tm, estimate]
    expect_equal(v, rep(unname(ols[tm]), length(v)), tolerance = 1e-8)
  }
})

test_that("FE mode recovers the within slope and tags gwfe + NFE", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwr_fe_panel()
  d$x_dm <- d$x - ave(d$x, d$unit); d$y_dm <- d$y - ave(d$y, d$unit)
  within <- unname(coef(lm(y_dm ~ x_dm, data = d))["x_dm"])

  out <- estimate_gwr(d, unit = "unit", formula = y ~ x,
    panel = "unit", time = "year", coords = c("lon", "lat"),
    kernel = "boxcar", adaptive = FALSE, bw = 1e6)

  expect_true(all(out$model_estimator == "gwfe"))
  xs <- out[term == "x", estimate]
  expect_equal(xs, rep(within, length(xs)), tolerance = 1e-6)
  expect_equal(unname(attr(out, "NFE")), length(unique(d$unit)))
})

test_that("fit_stats attaches finite within-R^2 and n_obs", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwr_fe_panel()
  out <- estimate_gwr(d, unit = "unit", formula = y ~ x,
    panel = "unit", time = "year", coords = c("lon", "lat"),
    kernel = "boxcar", adaptive = FALSE, bw = 1e6, fit_stats = TRUE, terms = "x")
  expect_true(all(c("r_squared", "n_obs") %in% names(out)))
  expect_true(all(is.finite(out$r_squared)))
  expect_true(all(out$r_squared >= 0 & out$r_squared <= 1))
})

test_that("standard errors are opt-in (off by default) and match lm under flat weights", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwr_fe_panel(); d$z <- rnorm(nrow(d))

  off <- estimate_gwr(d, unit = "unit", formula = y ~ x + z,
    coords = c("lon", "lat"), kernel = "boxcar", adaptive = FALSE, bw = 1e6)
  expect_true(all(is.na(off$se)))

  sm <- coef(summary(lm(y ~ x + z, data = d)))
  on <- estimate_gwr(d, unit = "unit", formula = y ~ x + z,
    coords = c("lon", "lat"), kernel = "boxcar", adaptive = FALSE, bw = 1e6,
    standard_errors = TRUE)
  one <- on[unit_id == "u1" & estimand == "mean"][match(rownames(sm), term)]
  expect_equal(one$se, unname(sm[, "Std. Error"]), tolerance = 1e-6)
})

test_that("terms filter restricts the returned terms", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .gwr_fe_panel(); d$z <- rnorm(nrow(d))
  out <- estimate_gwr(d, unit = "unit", formula = y ~ x + z,
    coords = c("lon", "lat"), kernel = "boxcar", adaptive = FALSE, bw = 1e6,
    terms = "x")
  expect_equal(unique(out$term), "x")
})
