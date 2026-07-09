# Shared, deterministic fixtures for the gwkit test suite.

# A regular n x n lattice of unit squares as an sf polygon layer (EPSG:4326).
# Centroids fall on (i - 0.5, j - 0.5), so the layer is also a clean Queen lattice.
make_sf_grid <- function(n = 5L, id_col = "pid") {
  skip_if_not_installed("sf")
  polys <- vector("list", n * n)
  ids   <- character(n * n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      m <- matrix(
        c(i - 1, j - 1,
          i,     j - 1,
          i,     j,
          i - 1, j,
          i - 1, j - 1),
        ncol = 2, byrow = TRUE
      )
      polys[[k]] <- sf::st_polygon(list(m))
      ids[k]     <- paste0("p", k)
      k <- k + 1L
    }
  }
  out <- sf::st_sf(id = ids, geometry = sf::st_sfc(polys, crs = 4326))
  names(out)[1] <- id_col
  out
}

# A p x p lattice of point units, each observed over `tt` with y = a + b*trend + noise.
make_point_panel <- function(p = 5L, tt = 0:4, seed = 1L) {
  set.seed(seed)
  grid <- expand.grid(ix = seq_len(p), iy = seq_len(p))
  n <- nrow(grid)
  b <- seq_len(n) / n
  do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      unit  = paste0("u", i),
      lon   = grid$ix[i],
      lat   = grid$iy[i],
      trend = tt,
      y     = 10 + b[i] * tt + stats::rnorm(length(tt), sd = 0.05),
      stringsAsFactors = FALSE
    )
  }))
}
