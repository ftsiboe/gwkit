# ============================================================
# Why does the panel fast path disagree with the direct path ONLY under
# Great Circle?
# ============================================================
# The fast path computes distances on UNIQUE locations and expands by index:
#
#     Du <- gw.dist(uxy,    tgt)      # n_uniq x n_tgt
#     De <- Du[uidx, ]                # n_obs  x n_tgt
#
# the direct path computes them on all observations:
#
#     Dd <- gw.dist(obs_xy, tgt)      # n_obs  x n_tgt
#
# These MUST be identical - a distance depends only on the two coordinates, and
# uxy[uidx, ] == obs_xy by construction. Euclidean, Manhattan and Minkowski
# agree to floating point. Great Circle does not, which means the assumption
# fails for longlat = TRUE and I do not know why.
#
# This script answers three questions, in order:
#   1. Is uxy[uidx, ] really equal to obs_xy?           (is the dedup right?)
#   2. Is De really equal to Dd?                        (is gw.dist pairwise?)
#   3. If not, WHERE do they differ, and by how much?
#
# Run from anywhere:
#   Rscript data-raw/scripts/diagnose_greatcircle_dedup.R
# ============================================================

suppressPackageStartupMessages({ library(GWmodel) })

set.seed(42)
n_loc <- 40; n_per <- 8
xy <- cbind(runif(n_loc, -100, -90), runif(n_loc, 35, 45))   # lon, lat
g  <- rep(seq_len(n_loc), each = n_per)
obs_xy <- xy[g, , drop = FALSE]
tgt_xy <- xy

# --- the fast path's dedup, verbatim from .gw_local_fit() -------------------
ukey   <- paste(obs_xy[, 1L], obs_xy[, 2L])
uidx   <- match(ukey, unique(ukey))
first  <- !duplicated(ukey)
uxy    <- obs_xy[first, , drop = FALSE]

cat("=== 1. is the dedup exact? ===\n")
cat("   n_obs =", nrow(obs_xy), " n_uniq =", nrow(uxy), "\n")
cat("   max |uxy[uidx, ] - obs_xy| =",
    max(abs(uxy[uidx, , drop = FALSE] - obs_xy)), "\n")
cat("   identical()                =",
    identical(uxy[uidx, , drop = FALSE], obs_xy), "\n\n")

cat("=== 2. does gw.dist agree between the two constructions? ===\n")
for (nm in c("Euclidean", "Great Circle")) {
  m  <- gwkit::resolve_distance_metric(nm)
  Du <- GWmodel::gw.dist(dp.locat = uxy,    rp.locat = tgt_xy, focus = 0,
                         p = m$p, theta = m$theta, longlat = m$longlat)
  Dd <- GWmodel::gw.dist(dp.locat = obs_xy, rp.locat = tgt_xy, focus = 0,
                         p = m$p, theta = m$theta, longlat = m$longlat)
  De <- Du[uidx, , drop = FALSE]

  cat(sprintf("\n   -- %s (longlat = %s, p = %s, theta = %s)\n",
              nm, m$longlat, m$p, m$theta))
  cat("      dim(Du) =", paste(dim(Du), collapse = "x"),
      "| dim(Dd) =", paste(dim(Dd), collapse = "x"),
      "| dim(De) =", paste(dim(De), collapse = "x"), "\n")
  cat("      class(Du) =", class(Du)[1], "| is.matrix =", is.matrix(Du), "\n")
  if (all(dim(De) == dim(Dd))) {
    d <- abs(De - Dd)
    cat("      max |De - Dd| =", max(d), "\n")
    cat("      # cells differing by > 1e-9:", sum(d > 1e-9), "/", length(d), "\n")
    if (max(d) > 1e-9) {
      w <- which(d == max(d), arr.ind = TRUE)[1, ]
      cat("      worst cell: obs row", w[1], "target", w[2], "\n")
      cat("        De =", De[w[1], w[2]], " Dd =", Dd[w[1], w[2]], "\n")
      cat("        obs coord :", obs_xy[w[1], ], "\n")
      cat("        uxy coord :", uxy[uidx[w[1]], ], "\n")
      cat("        tgt coord :", tgt_xy[w[2], ], "\n")
      # Is Dd self-consistent? Rows of the SAME location must be identical.
      j <- uidx[w[1]]
      rows <- which(uidx == j)
      cat("        rows of this location in Dd, target", w[2], ":\n          ",
          paste(signif(Dd[rows, w[2]], 10), collapse = " "), "\n")
      cat("        ^ if these are not all equal, gw.dist is not a pure pairwise\n",
          "          function of the two coordinates under longlat.\n")
    }
  } else {
    cat("      DIMENSIONS DIFFER - gw.dist is not returning dp x rp here.\n")
  }
}

cat("\n=== 3. does gw.weight see the same thing? ===\n")
m  <- gwkit::resolve_distance_metric("Great Circle")
Du <- GWmodel::gw.dist(uxy, tgt_xy, focus = 0, p = m$p, theta = m$theta, longlat = m$longlat)
Dd <- GWmodel::gw.dist(obs_xy, tgt_xy, focus = 0, p = m$p, theta = m$theta, longlat = m$longlat)
We <- matrix(GWmodel::gw.weight(Du[uidx, , drop = FALSE], bw = 120, kernel = "bisquare",
                                adaptive = TRUE), nrow = nrow(obs_xy))
Wd <- matrix(GWmodel::gw.weight(Dd, bw = 120, kernel = "bisquare",
                                adaptive = TRUE), nrow = nrow(obs_xy))
cat("   max |W_expanded - W_direct| =", max(abs(We - Wd)), "\n")
cat("   nonzero weights: expanded =", sum(We > 0), " direct =", sum(Wd > 0), "\n")
cat("\n   If the distances match but the weights do not, the problem is in\n")
cat("   gw.weight; if the distances already differ, it is in gw.dist.\n")
