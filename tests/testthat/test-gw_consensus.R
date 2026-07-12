# Deterministic lattice: the Queen-contiguity logic is exercised on a regular
# grid built in-test (GWmodel's datasets are not gridded, so a synthetic lattice
# is the appropriate fixture here).

make_lattice_class_dt <- function(nx = 4L, ny = 4L, settings = 3L,
                                   levels = c("A", "B")) {
  grid  <- expand.grid(ix = seq_len(nx), iy = seq_len(ny))
  units <- paste0("g", seq_len(nrow(grid)))
  truth <- rep(levels, length.out = nrow(grid))
  rows  <- lapply(seq_len(nrow(grid)), function(i) {
    data.frame(
      unit  = units[i],
      lon   = grid$ix[i],
      lat   = grid$iy[i],
      setting = paste0("s", seq_len(settings)),
      class = rep(truth[i], settings),      # all settings agree -> modal = truth
      stringsAsFactors = FALSE
    )
  })
  list(dt = do.call(rbind, rows), units = units, truth = truth)
}

test_that("gw_consensus_class (point) returns one row per unit with expected columns", {
  fx  <- make_lattice_class_dt()
  res <- gw_consensus_class(
    fx$dt, unit_col = "unit", class_col = "class",
    coords = c("lon", "lat"), class_levels = c("A", "B")
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), length(fx$units))
  expect_true(all(c("unit", "longitude", "latitude", "n_settings",
                    "modal_class", "modal_agreement",
                    "queen_class", "queen_order", "queen_agreement") %in% names(res)))
})

test_that("modal class is the mode across settings and agreement is 1 when unanimous", {
  fx  <- make_lattice_class_dt()
  res <- gw_consensus_class(
    fx$dt, unit_col = "unit", class_col = "class",
    coords = c("lon", "lat"), class_levels = c("A", "B")
  )
  expected <- fx$truth[match(res$unit, fx$units)]
  expect_identical(res$modal_class, expected)
  expect_true(all(res$modal_agreement == 1))
  expect_equal(res$n_settings, rep(3, nrow(res)))
})

test_that("queen vote resolves and reports the order it resolved at", {
  fx  <- make_lattice_class_dt()
  res <- gw_consensus_class(
    fx$dt, unit_col = "unit", class_col = "class",
    coords = c("lon", "lat"), class_levels = c("A", "B")
  )
  expect_true(all(!is.na(res$queen_class)))
  expect_true(all(res$queen_class %in% c("A", "B")))
  expect_true(all(res$queen_order >= 1L))
})
