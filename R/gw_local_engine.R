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
.gw_local_fit <- function(mm, y, obs_xy, tgt_xy,
                          p, theta, longlat, kernel, adaptive,
                          bw = NULL, bw_approach = "CV",
                          standard_errors = FALSE, fit_stats = FALSE,
                          chunk_size = 500L) {

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("Package 'GWmodel' is required.")

  mm <- as.matrix(mm); y <- as.numeric(y)
  k  <- ncol(mm); n_tgt <- nrow(tgt_xy)

  if (is.null(bw)) {
    d_self <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = obs_xy,
                               focus = 0, p = p, theta = theta, longlat = longlat)
    sub <- sp::SpatialPointsDataFrame(coords = obs_xy, data = data.frame(.y = y))
    bw  <- GWmodel::bw.gwr(.y ~ 1, data = sub, approach = bw_approach,
                           kernel = kernel, adaptive = adaptive,
                           p = p, theta = theta, longlat = longlat, dMat = d_self)
    rm(d_self)
  }

  coef <- matrix(NA_real_, n_tgt, k, dimnames = list(NULL, colnames(mm)))
  se   <- if (standard_errors) matrix(NA_real_, n_tgt, k) else NULL
  r2   <- if (fit_stats)       rep(NA_real_, n_tgt)       else NULL
  nobs <- rep(NA_integer_, n_tgt)

  chunks <- split(seq_len(n_tgt), ceiling(seq_len(n_tgt) / chunk_size))
  for (ch in chunks) {
    D <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy[ch, , drop = FALSE],
                          focus = 0, p = p, theta = theta, longlat = longlat)
    W <- GWmodel::gw.weight(vdist = D, bw = bw, kernel = kernel, adaptive = adaptive)
    W <- matrix(W, nrow = nrow(obs_xy), ncol = length(ch))

    for (jj in seq_along(ch)) {
      j <- ch[jj]; w <- W[, jj]
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
