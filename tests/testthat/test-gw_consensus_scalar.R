# Tests for the scalar (continuous) consensus in gw_consensus_scalar.R

# A tiny stacked table: 3 units (A, B, C) on a 1x3 lattice, 3 settings each.
#   A: 1, 2, 3   -> median 2,  all positive
#   B: -2,-4,-6  -> median -4, all negative
#   C: -1, 2, 3  -> median 2,  2/3 positive (sign disagreement)
.scalar_fixture <- function() {
  data.frame(
    grid_id   = rep(c("A", "B", "C"), each = 3),
    longitude = rep(c(1, 2, 3), each = 3),
    latitude  = rep(c(1, 1, 1), each = 3),
    kernel    = rep(c("k1", "k2", "k3"), times = 3),
    estimate  = c(1, 2, 3,  -2, -4, -6,  -1, 2, 3),
    stringsAsFactors = FALSE
  )
}

test_that("by_point returns the median consensus, sign agreement, and n_settings", {
  out <- gw_optimal_scalar_by_point(.scalar_fixture(), unit_col = "grid_id")
  out <- out[order(out$grid_id)]

  expect_s3_class(out, "data.table")
  expect_equal(out$grid_id, c("A", "B", "C"))
  expect_equal(out$n_settings, c(3L, 3L, 3L))
  expect_equal(out$consensus, c(2, -4, 2))          # per-unit medians
  expect_equal(out$sign_agreement, c(1, 1, round(2 / 3, 3)))
  expect_true(all(out$consensus_sd  >= 0))
  expect_true(all(out$consensus_lo <= out$consensus_hi))
  # expected column set, and no Queen columns unless requested
  expect_true(all(c("grid_id", "longitude", "latitude", "n_settings", "consensus",
                    "consensus_sd", "consensus_mad", "consensus_lo", "consensus_hi",
                    "sign_agreement") %in% names(out)))
  expect_false(any(c("queen_value", "queen_order", "queen_agreement") %in% names(out)))
})

test_that("agg_fun is honoured (mean, and a custom trimmed-mean closure)", {
  out_mean <- gw_optimal_scalar_by_point(.scalar_fixture(), "grid_id", agg_fun = mean)
  out_mean <- out_mean[order(out_mean$grid_id)]
  expect_equal(out_mean$consensus, c(2, -4, 4 / 3))          # per-unit means

  out_trim <- gw_optimal_scalar_by_point(.scalar_fixture(), "grid_id",
                                         agg_fun = function(x) mean(x, trim = 0.5))
  out_trim <- out_trim[order(out_trim$grid_id)]
  expect_equal(out_trim$consensus, c(2, -4, 2))              # 50% trim == median
})

test_that("probs control the reported spread quantiles", {
  out <- gw_optimal_scalar_by_point(.scalar_fixture(), "grid_id", probs = c(0, 1))
  out <- out[order(out$grid_id)]
  expect_equal(out$consensus_lo, c(1, -6, -1))              # per-unit minima
  expect_equal(out$consensus_hi, c(3, -2,  3))              # per-unit maxima
})

test_that("queen_smooth adds finite neighbourhood columns", {
  out <- gw_optimal_scalar_by_point(.scalar_fixture(), "grid_id", queen_smooth = TRUE)
  expect_true(all(c("queen_value", "queen_order", "queen_agreement") %in% names(out)))
  expect_true(all(is.finite(out$queen_value)))
  expect_true(all(out$queen_order >= 1L))
})

test_that("non-finite values are dropped before the summary", {
  df <- .scalar_fixture()
  df$estimate[1] <- NA_real_                                # A's first setting -> NA
  out <- gw_optimal_scalar_by_point(df, "grid_id")
  out <- out[order(out$grid_id)]
  expect_equal(out$n_settings[out$grid_id == "A"], 2L)
  expect_equal(out$consensus[out$grid_id == "A"], 2.5)     # median of c(2, 3)
})

test_that("by_polygon consensus over an sf grid, with correct shape and centroids", {
  skip_if_not_installed("sf")
  g   <- make_sf_grid(n = 3L, id_col = "pid")              # 9 polygons p1..p9
  ids <- g$pid
  # per polygon: 3 settings centred on v = index - 5 (spans negative..positive)
  df <- do.call(rbind, lapply(seq_along(ids), function(i) {
    v <- i - 5
    data.frame(pid = ids[i], kernel = c("k1", "k2", "k3"),
               estimate = c(v - 0.5, v, v + 0.5), stringsAsFactors = FALSE)
  }))

  out <- gw_optimal_scalar_by_polygon(df, unit_col = "pid", polygons = g,
                                      value_col = "estimate")
  out <- out[match(ids, out$pid)]

  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), length(ids))
  expect_equal(out$consensus, as.numeric(seq_along(ids) - 5))   # per-unit medians
  expect_true(all(is.finite(out$longitude) & is.finite(out$latitude)))

  out_q <- gw_optimal_scalar_by_polygon(df, "pid", polygons = g,
                                        value_col = "estimate", queen_smooth = TRUE)
  expect_true(all(c("queen_value", "queen_order", "queen_agreement") %in% names(out_q)))
  expect_true(all(is.finite(out_q$queen_value)))
})

test_that("by_polygon errors on a non-polygon input", {
  df <- data.frame(pid = "p1", kernel = "k1", estimate = 1)
  expect_error(
    gw_optimal_scalar_by_polygon(df, "pid", polygons = "not-a-polygon",
                                 value_col = "estimate")
  )
})
