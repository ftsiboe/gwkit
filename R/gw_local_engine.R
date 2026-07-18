# ============================================================
# Shared geographically weighted local-regression engine
# ============================================================
# One codepath for every gwkit local regression: distances from
# GWmodel::gw.dist(), kernel weights from GWmodel::gw.weight() (the SINGLE weight
# source across the package), and a per-focal weighted least-squares solve via
# stats::lm.wfit() - which is exactly what stats::lm() calls internally, so the
# point estimates are identical to a per-focal lm() but without the formula and
# summary overhead.
#
# By design this returns COEFFICIENTS ONLY. Standard errors are opt-in
# (`standard_errors = TRUE`); the intended inference path is the study's
# bootstrap, so SEs are off by default. A weighted within-R^2 (`fit_stats`) is
# available for model selection (e.g. per-county knot choice).
# ============================================================

# obs_xy, tgt_xy : n_obs x 2 and n_tgt x 2 coordinate matrices
# mm             : n_obs x k model matrix (built once, incl. intercept)
# y              : length-n_obs response
# returns: list(coef [n_tgt x k], se [n_tgt x k or NULL], r_squared [or NULL],
#               n_obs [per focal], bw)
# ============================================================
# Internal: gw.dist() with a GUARANTEED dp x rp orientation
# ============================================================
# GWmodel::gw.dist() is INCONSISTENT about the shape it returns:
#
#     longlat = FALSE (Minkowski path) -> dp x rp      (n_obs x n_targets)
#     longlat = TRUE  (great circle)   -> rp x dp      (n_targets x n_obs)   <-- !
#
# Measured, not inferred (data-raw/scripts/diagnose_greatcircle_dedup.R):
#
#     Euclidean    : dim(gw.dist(obs[320], tgt[40])) = 320 x 40
#     Great Circle : dim(gw.dist(obs[320], tgt[40])) =  40 x 320
#
# This silently corrupted every Great Circle fit. gw.weight() applies an adaptive
# bandwidth COLUMN-WISE, so with the transposed matrix each column was one
# observation's distances to the targets rather than one target's distances to
# the observations - the neighbourhood was taken along the wrong axis entirely.
# With bw = 120 and only 40 targets per column, every weight came back nonzero:
# 12,800 of 12,800, when a compact bisquare kernel should have produced ~4,480.
# The subsequent matrix(W, nrow = n_obs) then reshaped the transposed result and
# scrambled the observation-to-target mapping.
#
# Anything gwkit has produced under distance_metric = "Great Circle" is invalid.
# Other metrics are unaffected: they were dp x rp all along.
#
# Orientation is resolved by DIMENSION where that is unambiguous. When the matrix
# is square (n_obs == n_targets - the usual case here, since targets are the
# unique units) a second one-target call settles it. Note that when the targets
# ARE the observation locations the matrix is symmetric and orientation cannot
# matter; the probe only costs anything in that already-cheap case.
# ------------------------------------------------------------
.gw_dist_oriented <- function(dp, rp, p, theta, longlat) {
  n_dp <- nrow(dp); n_rp <- nrow(rp)
  D <- GWmodel::gw.dist(dp.locat = dp, rp.locat = rp, focus = 0,
                        p = p, theta = theta, longlat = longlat)
  D <- as.matrix(D)

  if (nrow(D) == n_dp && ncol(D) == n_rp) {
    if (n_dp != n_rp) return(D)                      # unambiguous: already dp x rp
    # Square: ask for a single target and see which way it comes back.
    pr <- as.matrix(GWmodel::gw.dist(dp.locat = dp, rp.locat = rp[1L, , drop = FALSE],
                                     focus = 0, p = p, theta = theta, longlat = longlat))
    if (n_dp > 1L && nrow(pr) == 1L && ncol(pr) == n_dp) return(t(D))
    return(D)
  }
  if (nrow(D) == n_rp && ncol(D) == n_dp) return(t(D))

  stop("GWmodel::gw.dist() returned an unexpected shape: ",
       paste(dim(D), collapse = " x "), " for dp = ", n_dp, ", rp = ", n_rp)
}


# ============================================================
# Internal: per-location sufficient statistics (PANEL fast path)
# ============================================================
# When several observations share a location - a panel, where each unit appears
# once per period at IDENTICAL coordinates - the geographic weight is constant
# within a unit. The weighted normal equations then collapse from a sum over
# OBSERVATIONS to a sum over LOCATIONS:
#
#   X'WX = sum_i w_i x_i x_i'                 (n_obs terms)
#        = sum_j w_j sum_{t in j} x_jt x_jt'  (w constant within unit j)
#        = sum_j w_j S_j                      (n_uniq terms)
#
# and likewise X'Wy = sum_j w_j s_j. Everything .gw_local_fit() returns follows
# from the same per-unit quantities:
#
#   rss  = sum_j w_j [ q_j - 2 b's_j + b'S_j b ]
#   tss  = sum_j w_j [ q_j - 2 ybar m_j + n_j ybar^2 ],  ybar = sum_j w_j m_j / sum_j w_j n_j
#   nobs = sum_j n_j [w_j > 0]
#
# This is an ALGEBRAIC IDENTITY, not an approximation: verified to ~2e-14
# relative error against the direct lm.wfit path over random weightings,
# including unbalanced panels and weights that are zero outside a neighbourhood.
#
# Why it matters: the direct path calls lm.wfit() once per focal point, each
# factorising the FULL n_obs x k model matrix. For a 13,626-grid x 79-year panel
# that is 1.08M rows re-factorised 13,626 times per model - measured at ~6.5h for
# an 8-year window and extrapolating to ~64h at 79 years. The fast path replaces
# every focal's solve with a slice of one chunk-wide matmul over 13,626 rows.
# ------------------------------------------------------------
.gw_panel_stats <- function(mm, y, uidx, n_uniq) {
  k   <- ncol(mm)
  one <- rep(1, length(y))
  # rowsum() groups by uidx, whose values are 1..n_uniq, so rows come back in
  # unit order and align with uxy.
  S <- matrix(0, n_uniq, k * k)
  for (a in seq_len(k)) for (b in seq_len(k))
    S[, (a - 1L) * k + b] <- rowsum(mm[, a] * mm[, b], uidx, reorder = TRUE)[, 1L]
  sv <- matrix(0, n_uniq, k)
  for (a in seq_len(k))
    sv[, a] <- rowsum(mm[, a] * y, uidx, reorder = TRUE)[, 1L]
  list(S = S, s = sv, k = k,
       q = rowsum(y * y, uidx, reorder = TRUE)[, 1L],
       m = rowsum(y,     uidx, reorder = TRUE)[, 1L],
       n = rowsum(one,   uidx, reorder = TRUE)[, 1L])
}


# obs_xy, tgt_xy : n_obs x 2 and n_tgt x 2 coordinate matrices
# mm             : n_obs x k model matrix (built once, incl. intercept)
# y              : length-n_obs response
# returns: list(coef [n_tgt x k], se [n_tgt x k or NULL], r_squared [or NULL],
#               n_obs [per focal], bw)
# `.force_direct` exists only so the equivalence test can drive BOTH paths on the
# same inputs; the panel path is an algebraic identity, and the test is what
# demonstrates that rather than asserts it. Never set it in production - it is
# the slow path.
.gw_local_fit <- function(mm, y, obs_xy, tgt_xy,
                          p, theta, longlat, kernel, adaptive,
                          bw = NULL, bw_approach = "CV",
                          standard_errors = FALSE, fit_stats = FALSE,
                          chunk_size = 500L, .force_direct = FALSE) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required.")

  mm <- as.matrix(mm); y <- as.numeric(y)
  k  <- ncol(mm); n_tgt <- nrow(tgt_xy); n_obs <- nrow(obs_xy)

  if (is.null(bw)) {
    # Self-distances: dp == rp, so this matrix is square and symmetric and the
    # gw.dist() orientation bug (see .gw_dist_oriented) cannot bite here.
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(coords = obs_xy, data = data.frame(.y = y))
    bw  <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                           kernel = kernel, adaptive = adaptive,
                           p = p, theta = theta, longlat = longlat, dMat = d_self)
    rm(d_self)
  }

  # --- detect repeated locations -------------------------------------------
  ukey   <- paste(obs_xy[, 1L], obs_xy[, 2L])
  uidx   <- match(ukey, unique(ukey))          # obs -> unit id (1..n_uniq)
  first  <- !duplicated(ukey)
  uxy    <- obs_xy[first, , drop = FALSE]
  n_uniq <- nrow(uxy)
  panel  <- (n_uniq < n_obs) && !isTRUE(.force_direct)   # nothing to gain when all unique

  st <- if (panel) .gw_panel_stats(mm, y, uidx, n_uniq) else NULL

  coef <- matrix(NA_real_, n_tgt, k, dimnames = list(NULL, colnames(mm)))
  se   <- if (standard_errors) matrix(NA_real_, n_tgt, k) else NULL
  r2   <- if (fit_stats)       rep(NA_real_, n_tgt)       else NULL
  nobs <- rep(NA_integer_, n_tgt)

  # Distances are computed on UNIQUE locations only - in a panel the period
  # replicates sit on top of each other, so the direct path recomputes every
  # distance n_periods times over.
  dist_xy <- if (panel) uxy else obs_xy

  # Cap the transient D/W allocation. The direct path fixed chunk_size at 500,
  # which at n_obs = 1.08M is 2 x 4.3 GB per worker - times 228 workers that is
  # 1.79 TB and the OOM killer takes the run apart. Size the chunk to the data
  # instead of the other way round.
  eff_chunk <- if (panel) {
    max(1L, min(chunk_size, floor(256 * 1024^2 / (8 * n_obs))))
  } else chunk_size

  chunks <- split(seq_len(n_tgt), ceiling(seq_len(n_tgt) / eff_chunk))
  for (ch in chunks) {
    # .gw_dist_oriented(), not gw.dist(): the raw call returns rp x dp under
    # longlat = TRUE. See the note on that function.
    Du <- .gw_dist_oriented(dist_xy, tgt_xy[ch, , drop = FALSE], p, theta, longlat)

    if (!panel) {
      W <- GWmodel::gw.weight(vdist = Du, bw = bw, kernel = kernel, adaptive = adaptive)
      W <- matrix(W, nrow = n_obs, ncol = length(ch))
    } else {
      # Expand to observations ONLY to compute the weights, then collapse back.
      # This keeps gw.weight()'s semantics byte-for-byte identical to the direct
      # path - which matters for `adaptive`, where the bandwidth is the bw-th
      # smallest distance among OBSERVATIONS, not among locations. Deriving that
      # threshold from unique distances plus multiplicities would be faster
      # still, but would re-implement GWmodel's convention rather than call it.
      De <- Du[uidx, , drop = FALSE]
      We <- GWmodel::gw.weight(vdist = De, bw = bw, kernel = kernel, adaptive = adaptive)
      We <- matrix(We, nrow = n_obs, ncol = length(ch))
      W  <- We[first, , drop = FALSE]            # n_uniq x |ch|: one weight per unit
      rm(De, We)
    }

    if (panel) {
      # One matmul per chunk replaces |ch| calls to lm.wfit().
      XtWX_all <- crossprod(W, st$S)             # |ch| x k^2
      XtWy_all <- crossprod(W, st$s)             # |ch| x k
      qw <- as.vector(crossprod(W, st$q))
      mw <- as.vector(crossprod(W, st$m))
      nw <- as.vector(crossprod(W, st$n))        # sum of weights over OBSERVATIONS
    }

    for (jj in seq_along(ch)) {
      j <- ch[jj]; w <- W[, jj]

      if (!panel) {
        nobs[j] <- sum(w > 0)
        fit <- tryCatch(stats::lm.wfit(x = mm, y = y, w = w), error = function(e) NULL)
        if (is.null(fit)) next
        coef[j, ] <- fit$coefficients

        if (fit_stats || standard_errors) {
          r   <- fit$residuals
          rss <- sum(w * r^2)
          if (fit_stats) {
            ybar <- sum(w * y) / sum(w)
            tss  <- sum(w * (y - ybar)^2)
            r2[j] <- if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
          }
          if (standard_errors) {
            dfres <- sum(w > 0) - fit$rank
            if (is.finite(dfres) && dfres > 0) {
              tryCatch({
                XtWX <- crossprod(mm, mm * w)             # X'WX (k x k)
                inv  <- solve(XtWX)                        # (X'WX)^-1
                se[j, ] <- sqrt((rss / dfres) * diag(inv)) # Var(b)=sigma^2 (X'WX)^-1
              }, error = function(e) {})
            }
          }
        }

      } else {
        # byrow: column (a-1)*k + b of S holds sum_t x_a x_b.
        XtWX <- matrix(XtWX_all[jj, ], k, k, byrow = TRUE)
        XtWy <- XtWy_all[jj, ]
        nobs[j] <- as.integer(sum(st$n[w > 0]))            # in OBSERVATIONS, as before

        # ------------------------------------------------------------
        # RANK CHECK - do not skip this.
        # ------------------------------------------------------------
        # lm.wfit() (the direct path) solves via a PIVOTED QR with a rank
        # tolerance: when a local neighbourhood cannot identify a coefficient it
        # drops the column and returns NA. solve() has no such check and will
        # return plausible-looking numbers for a near-singular system.
        #
        # This is not hypothetical. The equivalence test caught exactly this
        # under the Great Circle metric: the direct path reported
        # coef = (3.09, NA) with r^2 = 0 - trend not identified - while this path
        # happily returned (2.55, 0.49). The fabricated value was even close to
        # the simulation truth, which is precisely why it is dangerous: those
        # grids would have entered the trend classification as real estimates
        # instead of NA.
        #
        # rcond() on the k x k normal matrix is the cheap guard. Note X'WX squares
        # the conditioning of X, so compare against sqrt(eps) rather than eps.
        # Leaving the row as NA is conservative: it refuses where lm.wfit would
        # have returned a partially-NA row, so this path never claims MORE than
        # the direct path does.
        # ------------------------------------------------------------
        rc <- tryCatch(rcond(XtWX), error = function(e) 0)
        if (!is.finite(rc) || rc < sqrt(.Machine$double.eps)) next

        b <- tryCatch(solve(XtWX, XtWy), error = function(e) NULL)
        if (is.null(b) || !all(is.finite(b))) next
        coef[j, ] <- b

        if (fit_stats || standard_errors) {
          rss <- qw[jj] - 2 * sum(b * XtWy) + as.vector(t(b) %*% XtWX %*% b)
          if (fit_stats) {
            ybar  <- mw[jj] / nw[jj]
            tss   <- qw[jj] - 2 * ybar * mw[jj] + ybar^2 * nw[jj]
            r2[j] <- if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
          }
          if (standard_errors) {
            dfres <- nobs[j] - k
            if (is.finite(dfres) && dfres > 0) {
              tryCatch({
                inv <- solve(XtWX)
                se[j, ] <- sqrt((rss / dfres) * diag(inv))
              }, error = function(e) {})
            }
          }
        }
      }
    }
  }
  list(coef = coef, se = se, r_squared = r2, n_obs = nobs, bw = bw)
}


# ============================================================
# Multi-kernel companion of .gw_local_fit()
# ============================================================
# Fits the SAME model matrix under several kernels while computing the
# (kernel-independent) obs->target distance matrix only ONCE per target chunk and
# reusing it across kernels - the regression analogue of estimate_gwss_kernels()'s
# distance reuse. Memory profile is identical to .gw_local_fit() (distances are
# still chunked over targets, never materialised whole), so this is safe at full
# scale; the saving is the eliminated redundant gw.dist() calls.
#
# `bw` may be NULL (select once), a single numeric (used for every kernel), or a
# named numeric keyed by kernel (per-kernel reuse). When NULL, `bandwidth`
# controls selection: "shared" runs ONE bw.gwr() (on kernels[1]) and reuses it
# for all kernels; "per_kernel" runs bw.gwr() once per kernel. Either way the
# self-distance matrix used for selection is built once.
#
# Returns list(fits = named list per kernel of
#   list(coef, se, r_squared, n_obs, bw), bw = named numeric of bandwidths).
.gw_local_fit_kernels <- function(mm, y, obs_xy, tgt_xy,
                                  p, theta, longlat, kernels, adaptive,
                                  bw = NULL, bandwidth = "shared", bw_approach = "CV",
                                  standard_errors = FALSE, fit_stats = FALSE,
                                  chunk_size = 500L) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required.")

  mm <- as.matrix(mm); y <- as.numeric(y)
  k  <- ncol(mm); n_tgt <- nrow(tgt_xy); nk <- length(kernels)

  # --- bandwidth(s), computed once, keyed by kernel --------------------------
  if (is.null(bw)) {
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(coords = obs_xy, data = data.frame(.y = y))
    if (identical(bandwidth, "shared")) {
      b0 <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                            kernel = kernels[[1L]], adaptive = adaptive,
                            p = p, theta = theta, longlat = longlat, dMat = d_self)
      bw_named <- stats::setNames(rep(b0, nk), kernels)
    } else {
      bw_named <- stats::setNames(vapply(kernels, function(km)
        GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                        kernel = km, adaptive = adaptive,
                        p = p, theta = theta, longlat = longlat, dMat = d_self),
        numeric(1)), kernels)
    }
    rm(d_self)
  } else if (is.null(names(bw))) {
    bw_named <- stats::setNames(rep(as.numeric(bw), length.out = nk), kernels)
  } else {
    bw_named <- bw[kernels]
    if (any(is.na(bw_named)))
      stop("named `bw` is missing entries for kernel(s): ",
           paste(kernels[is.na(bw_named)], collapse = ", "))
  }

  # --- per-kernel output stores ----------------------------------------------
  mk    <- function() matrix(NA_real_, n_tgt, k, dimnames = list(NULL, colnames(mm)))
  coefL <- stats::setNames(lapply(kernels, function(i) mk()), kernels)
  seL   <- if (standard_errors) stats::setNames(lapply(kernels, function(i) mk()), kernels) else NULL
  r2L   <- if (fit_stats)       stats::setNames(lapply(kernels, function(i) rep(NA_real_, n_tgt)), kernels) else NULL
  nobsL <- stats::setNames(lapply(kernels, function(i) rep(NA_integer_, n_tgt)), kernels)

  chunks <- split(seq_len(n_tgt), ceiling(seq_len(n_tgt) / chunk_size))
  for (ch in chunks) {
    D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy[ch, , drop = FALSE],
                          focus = 0, p = p, theta = theta, longlat = longlat)   # once per chunk
    for (km in kernels) {
      W <- GWmodel::gw.weight(vdist = D, bw = bw_named[[km]], kernel = km, adaptive = adaptive)
      W <- matrix(W, nrow = nrow(obs_xy), ncol = length(ch))
      for (jj in seq_along(ch)) {
        j <- ch[jj]; w <- W[, jj]
        nobsL[[km]][j] <- sum(w > 0)
        fit <- tryCatch(stats::lm.wfit(x = mm, y = y, w = w), error = function(e) NULL)
        if (is.null(fit)) next
        coefL[[km]][j, ] <- fit$coefficients

        if (fit_stats || standard_errors) {
          r   <- fit$residuals
          rss <- sum(w * r^2)
          if (fit_stats) {
            ybar <- sum(w * y) / sum(w)
            tss  <- sum(w * (y - ybar)^2)
            r2L[[km]][j] <- if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
          }
          if (standard_errors) {
            dfres <- sum(w > 0) - fit$rank
            if (is.finite(dfres) && dfres > 0) {
              tryCatch({
                XtWX <- crossprod(mm, mm * w)
                inv  <- solve(XtWX)
                seL[[km]][j, ] <- sqrt((rss / dfres) * diag(inv))
              }, error = function(e) {})
            }
          }
        }
      }
    }
  }

  fits <- stats::setNames(lapply(kernels, function(km) list(
    coef      = coefL[[km]],
    se        = if (!is.null(seL)) seL[[km]] else NULL,
    r_squared = if (!is.null(r2L)) r2L[[km]] else NULL,
    n_obs     = nobsL[[km]],
    bw        = bw_named[[km]])), kernels)
  list(fits = fits, bw = bw_named)
}
