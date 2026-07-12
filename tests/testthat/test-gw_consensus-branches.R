# Exercises the gw_consensus branches the lattice test does not reach:
# the polygon variant (shared-boundary Queen contiguity), include_self = FALSE,
# tie-driven Queen order expansion, and isolated units.

test_that("gw_consensus_class (polygon) works on an sf polygon layer", {
  skip_if_not_installed("sf")

  pg    <- make_sf_grid(4L, id_col = "pid")
  ids   <- pg$pid
  truth <- rep(c("A", "B"), length.out = length(ids))

  cls <- do.call(rbind, lapply(seq_along(ids), function(i) {
    data.frame(pid = ids[i], setting = paste0("s", 1:3),
               class = rep(truth[i], 3), stringsAsFactors = FALSE)
  }))

  res <- gw_consensus_class(
    cls, unit_col = "pid", geometry =pg, class_col = "class",
    poly_id = "pid", class_levels = c("A", "B")
  )

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), length(ids))
  expect_true(all(c("pid", "longitude", "latitude", "modal_class",
                    "queen_class", "queen_order") %in% names(res)))
  expect_identical(res$modal_class, truth[match(res$pid, ids)])
  expect_true(all(res$modal_agreement == 1))
})

test_that("gw_consensus_class (polygon) rejects a non-spatial polygons argument", {
  cls <- data.frame(pid = "p1", class = "A", stringsAsFactors = FALSE)
  expect_error(
    gw_consensus_class(cls, unit_col = "pid",
                                geometry =data.frame(pid = "p1"),
                                class_col = "class"),
    "SpatVector or an sf"
  )
})

test_that("Queen vote widens the order to break a tie, and isolates return NA", {
  n      <- 5L
  grid   <- expand.grid(ix = seq_len(n), iy = seq_len(n))
  units  <- paste0("g", seq_len(nrow(grid)))
  centre <- which(grid$ix == 3 & grid$iy == 3)

  # first-order Queen neighbours of the centre cell
  nbr <- which(abs(grid$ix - 3) <= 1 & abs(grid$iy - 3) <= 1 &
                 seq_len(nrow(grid)) != centre)

  cls          <- rep("A", nrow(grid))   # ring 2 is all "A"
  cls[nbr[5:8]] <- "B"                   # 4 A / 4 B among first-order neighbours
  cls[centre]  <- "B"                    # ignored: include_self = FALSE

  dt <- rbind(
    data.frame(unit = units, lon = grid$ix, lat = grid$iy, class = cls,
               stringsAsFactors = FALSE),
    # an isolated cell on the same lattice step, with no Queen neighbours
    data.frame(unit = "iso", lon = 10, lat = 10, class = "B",
               stringsAsFactors = FALSE)
  )

  res <- gw_consensus_class(
    dt, unit_col = "unit", class_col = "class", coords = c("lon", "lat"),
    class_levels = c("A", "B"), include_self = FALSE
  )

  # the centre ties at first order (4 vs 4) and only resolves after widening
  expect_gt(res$queen_order[res$unit == units[centre]], 1L)
  expect_identical(res$queen_class[res$unit == units[centre]], "A")

  # no neighbours and no self vote -> no Queen class
  expect_true(is.na(res$queen_class[res$unit == "iso"]))

  # the modal class is per-unit and unaffected by the Queen step
  expect_identical(res$modal_class[res$unit == "iso"], "B")
})
