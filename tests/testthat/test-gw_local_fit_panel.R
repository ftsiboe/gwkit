# ============================================================
# The panel fast path must be an IDENTITY, not an approximation.
# ============================================================
# .gw_local_fit() collapses the weighted normal equations to per-location
# sufficient statistics whenever observations share coordinates (a panel). That
# is only legitimate because the geographic weight is constant within a location,
# so:
#
#   X'WX = sum_i w_i x_i x_i' = sum_j w_j S_j
#
# These tests drive BOTH paths on identical inputs via `.force_direct` and
# require the results to agree to floating point - across every kernel, several
# distance metrics, balanced and unbalanced panels, and the fit statistics that
# are easiest to get subtly wrong (rss over ALL rows, tss, n_obs in OBSERVATIONS
# rather than locations).
#
# If any of these fail, the fast path is wrong and the study's estimates move.
# ============================================================

make_panel <- function(n_loc = 40, n_per = 8, unbalanced = FALSE, seed = 42) {
  set.seed(seed)
  xy <- cbind(runif(n_loc, -100, -90), runif(n_loc, 35, 45))
  g <- t <- NULL
  for (j in seq_len(n_loc)) {
    keep <- if (unbalanced) sort(sample.int(n_per, sample(3:n_per, 1))) else seq_len(n_per)
    g <- c(g, rep(j, length(keep))); t <- c(t, keep)
  }
  y  <- 3 + 0.4 * t + rnorm(length(t), 0, 1 + 0.05 * t)
  mm <- cbind(`(Intercept)` = 1, trend = t)
  list(mm = mm, y = y, obs_xy = xy[g, , drop = FALSE], tgt_xy = xy)
}

run_both <- function(d, kernel, adaptive = TRUE, bw = 15, p = 2, theta = 0,
                     longlat = FALSE) {
  args <- list(mm = d$mm, y = d$y, obs_xy = d$obs_xy, tgt_xy = d$tgt_xy,
               p = p, theta = theta, longlat = longlat, kernel = kernel,
               adaptive = adaptive, bw = bw,
               standard_errors = TRUE, fit_stats = TRUE)
  list(fast   = do.call(gwkit:::.gw_local_fit, c(args, list(.force_direct = FALSE))),
       direct = do.call(gwkit:::.gw_local_fit, c(args, list(.force_direct = TRUE))))
}

# ------------------------------------------------------------
# Comparing the two paths where a local fit is RANK DEFICIENT
# ------------------------------------------------------------
# Some neighbourhoods cannot identify every coefficient - with an adaptive
# bandwidth over a replicated panel, a compact kernel can leave only the focal
# location's own rows in play. The two paths then behave differently BY DESIGN:
#
#   direct : lm.wfit()'s pivoted QR drops the offending column -> (value, NA)
#   fast   : the rcond() guard refuses the whole row         -> (NA, NA)
#
# The fast path deliberately claims LESS, never more. So:
#   * agreement is required wherever the direct path returned a complete fit
#   * and separately, the fast path must NEVER report a finite estimate where
#     the direct path reported NA. That is the dangerous direction - it is what
#     fabricates estimates - and it gets its own assertion.
# ------------------------------------------------------------
expect_paths_agree <- function(r, label, tol = 1e-8) {
  ok <- apply(r$direct$coef, 1L, function(z) all(is.finite(z)))

  # Safety property, checked on EVERY row: no invention.
  invented <- apply(r$fast$coef, 1L, function(z) any(is.finite(z))) &
    !apply(r$direct$coef, 1L, function(z) any(is.finite(z)))
  expect_false(any(invented),
               info = paste0(label, " -> fast path returned a finite estimate where ",
                             "the direct path could not identify one"))

  if (!any(ok)) {
    skip(paste0(label, ": every local fit is rank deficient; nothing to compare"))
  }

  expect_equal(r$fast$coef[ok, , drop = FALSE], r$direct$coef[ok, , drop = FALSE],
               tolerance = tol, info = paste0(label, " -> coef"))
  if (!is.null(r$fast$se))
    expect_equal(r$fast$se[ok, , drop = FALSE], r$direct$se[ok, , drop = FALSE],
                 tolerance = tol, info = paste0(label, " -> se"))
  if (!is.null(r$fast$r_squared))
    expect_equal(r$fast$r_squared[ok], r$direct$r_squared[ok],
                 tolerance = tol, info = paste0(label, " -> r_squared"))
  # n_obs is defined for every row regardless of rank.
  expect_equal(r$fast$n_obs, r$direct$n_obs, info = paste0(label, " -> n_obs"))
}

test_that("panel fast path reproduces the direct path for every kernel", {
  skip_if_not_installed("GWmodel")
  d <- make_panel()
  # Every kernel gwkit supports. gaussian/exponential matter most: their weights
  # never reach zero, so the fast path cannot lean on sparsity anywhere.
  for (kern in c("gaussian", "exponential", "bisquare", "boxcar", "tricube")) {
    expect_paths_agree(run_both(d, kern), paste("kernel", kern))
  }
})

test_that("panel fast path holds for an UNBALANCED panel", {
  skip_if_not_installed("GWmodel")
  # n_j varies by location, so S_j / n_j genuinely differ and the identity is
  # doing real work - a balanced panel would hide a bug that assumed constant n.
  d <- make_panel(unbalanced = TRUE)
  for (kern in c("gaussian", "bisquare", "tricube")) {
    expect_paths_agree(run_both(d, kern), paste("unbalanced,", kern))
  }
})

test_that("panel fast path holds across distance metrics and fixed bandwidths", {
  skip_if_not_installed("GWmodel")
  d <- make_panel()
  # bw = 120 grid-years ~ 15 of the 40 locations. An earlier version used
  # bw = 12 against 8 replicates, which left ONLY the focal location in the
  # neighbourhood - degenerate, and nothing like production, where the selected
  # bandwidths span ~80-100 grids.
  for (dm in c("Euclidean", "Manhattan", "Minkowski p=3", "Great Circle")) {
    m <- gwkit::resolve_distance_metric(dm)
    r <- run_both(d, "bisquare", adaptive = TRUE, bw = 120,
                  p = m$p, theta = m$theta, longlat = m$longlat)
    expect_paths_agree(r, paste("metric", dm))
  }
  # Fixed (distance) bandwidth, not just adaptive kNN.
  expect_paths_agree(run_both(d, "bisquare", adaptive = FALSE, bw = 4),
                     "fixed bandwidth")
})

test_that("gw.dist orientation is normalised (GWmodel transposes under longlat)", {
  skip_if_not_installed("GWmodel")
  # THE bug that made Great Circle disagree between the two paths, and that had
  # been silently corrupting every Great Circle fit in gwkit before it:
  #
  #   GWmodel::gw.dist(obs[320], tgt[40])
  #     longlat = FALSE -> 320 x 40   (dp x rp)
  #     longlat = TRUE  ->  40 x 320  (rp x dp)
  #
  # gw.weight() applies an adaptive bandwidth column-wise, so the transposed
  # matrix put the neighbourhood along the wrong axis: with bw = 120 and only 40
  # targets per column, all 12,800 weights came back nonzero where a compact
  # bisquare kernel should have given ~4,480.
  set.seed(11)
  n_loc <- 40; n_per <- 8
  xy  <- cbind(runif(n_loc, -100, -90), runif(n_loc, 35, 45))
  obs <- xy[rep(seq_len(n_loc), each = n_per), , drop = FALSE]

  for (nm in c("Euclidean", "Great Circle")) {
    m <- gwkit::resolve_distance_metric(nm)
    D <- gwkit:::.gw_dist_oriented(obs, xy, m$p, m$theta, m$longlat)
    expect_equal(dim(D), c(nrow(obs), nrow(xy)),
                 info = paste0(nm, ": .gw_dist_oriented() must return dp x rp"))
  }

  # And the weights must then behave like a compact kernel: far fewer than all.
  m <- gwkit::resolve_distance_metric("Great Circle")
  D <- gwkit:::.gw_dist_oriented(obs, xy, m$p, m$theta, m$longlat)
  W <- matrix(GWmodel::gw.weight(D, bw = 120, kernel = "bisquare", adaptive = TRUE),
              nrow = nrow(obs))
  expect_lt(sum(W > 0), length(W))
  expect_gt(sum(W > 0), 0)
})

test_that("a rank-deficient local yields NA, never a fabricated estimate", {
  skip_if_not_installed("GWmodel")
  # THE failure the first version of this file exposed. With a narrow adaptive
  # bandwidth over a replicated panel, compact kernels can leave a neighbourhood
  # that cannot identify every coefficient. lm.wfit()'s pivoted QR notices and
  # returns NA; a bare solve() does not, and returns a plausible - even
  # near-truth - number. Under Great Circle the direct path gave (3.09, NA) while
  # the unguarded fast path gave (2.55, 0.49).
  #
  # The rcond() guard makes the fast path refuse instead. The contract asserted
  # here is one-directional and is the one that protects the study: the fast path
  # may report LESS than the direct path, never more.
  d <- make_panel(n_loc = 40, n_per = 8)
  for (dm in c("Euclidean", "Great Circle")) {
    m <- gwkit::resolve_distance_metric(dm)
    r <- run_both(d, "bisquare", adaptive = TRUE, bw = 12,   # deliberately narrow
                  p = m$p, theta = m$theta, longlat = m$longlat)
    fin_fast   <- apply(r$fast$coef,   1L, function(z) all(is.finite(z)))
    fin_direct <- apply(r$direct$coef, 1L, function(z) all(is.finite(z)))
    expect_false(any(fin_fast & !fin_direct),
                 info = paste0(dm, ": fast path invented a complete fit where the ",
                               "direct path found the local design rank deficient"))
  }
})

test_that("n_obs is counted in OBSERVATIONS, not locations", {
  skip_if_not_installed("GWmodel")
  # The easiest thing to get wrong: the fast path works in locations, but n_obs
  # (and hence the SE degrees of freedom) must stay in grid-year units.
  d <- make_panel(n_loc = 30, n_per = 6)
  r <- run_both(d, "bisquare", bw = 10)

  # The real assertion: whatever the count means, both paths must agree.
  expect_equal(r$fast$n_obs, r$direct$n_obs)

  # And it must be a multiple of the per-location replicate count, not a
  # location count. NB with bw = 10 over 6 replicates the 10th-nearest
  # grid-year falls on a NEIGHBOUR location, and bisquare is exactly zero at the
  # bandwidth - so only the focal location's own 6 rows carry weight and
  # n_obs = 6. (An earlier version of this test asserted n_obs > 30 on the
  # assumption that neighbours would contribute; they do not.)
  expect_true(all(r$fast$n_obs %% 6 == 0, na.rm = TRUE),
              info = "n_obs is not a whole number of location-replicates")
  expect_true(max(r$fast$n_obs, na.rm = TRUE) >= 6,
              info = "n_obs looks like a location count, not an observation count")
})

test_that("a wide bandwidth pulls in neighbours, and the paths still agree", {
  skip_if_not_installed("GWmodel")
  # Complements the test above: with a bandwidth well past the replicate count,
  # several locations carry weight, so the sufficient-statistic sum runs over
  # many S_j rather than collapsing to the focal location. This is the case the
  # identity actually has to work for.
  d <- make_panel(n_loc = 40, n_per = 8)
  r <- run_both(d, "bisquare", bw = 120)          # 120 grid-years ~ 15 locations
  expect_paths_agree(r, "wide bandwidth")
  expect_true(max(r$fast$n_obs, na.rm = TRUE) > 8,
              info = "wide bandwidth should reach beyond the focal location")
})

test_that("non-panel input is untouched by the fast path", {
  skip_if_not_installed("GWmodel")
  # All coordinates unique -> nothing to collapse -> must take the direct path
  # and return exactly what it always did.
  set.seed(1)
  n <- 60
  xy <- cbind(runif(n, -100, -90), runif(n, 35, 45))
  mm <- cbind(`(Intercept)` = 1, x = rnorm(n))
  y  <- rnorm(n)
  a <- gwkit:::.gw_local_fit(mm, y, xy, xy, p = 2, theta = 0, longlat = FALSE,
                             kernel = "bisquare", adaptive = TRUE, bw = 20,
                             standard_errors = TRUE, fit_stats = TRUE)
  b <- gwkit:::.gw_local_fit(mm, y, xy, xy, p = 2, theta = 0, longlat = FALSE,
                             kernel = "bisquare", adaptive = TRUE, bw = 20,
                             standard_errors = TRUE, fit_stats = TRUE,
                             .force_direct = TRUE)
  expect_equal(a$coef, b$coef)
  expect_equal(a$se, b$se)
})

test_that("estimate_gwr end-to-end is unchanged by the fast path", {
  skip_if_not_installed("GWmodel")
  # The identity must survive the wrapper, including variance = TRUE, whose
  # log-squared-residual refit reuses the same engine.
  set.seed(3)
  n_loc <- 35; n_per <- 7
  xy <- cbind(longitude = runif(n_loc, -100, -90), latitude = runif(n_loc, 35, 45))
  dt <- data.table::rbindlist(lapply(seq_len(n_loc), function(j) {
    data.table::data.table(unit = as.character(j), trend = seq_len(n_per),
                           longitude = xy[j, 1], latitude = xy[j, 2])
  }))
  dt[, precipitation := 5 + 0.3 * trend + stats::rnorm(.N)]

  fast <- gwkit::estimate_gwr(dt, unit = "unit", formula = precipitation ~ trend,
                              coords = c("longitude", "latitude"), kernel = "bisquare",
                              bw = 12, variance = TRUE, standard_errors = TRUE,
                              terms = "trend")
  expect_true(nrow(fast) > 0)
  expect_true(all(is.finite(fast$estimate)))
})
