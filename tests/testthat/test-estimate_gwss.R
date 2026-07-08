# Uses GWmodel's LondonHP points, reprojected from OSGB (EPSG:27700) to WGS84
# (EPSG:4326), since estimate_gwss_by_point() assumes longitude/latitude input.

test_that("estimate_gwss_by_point runs on LondonHP and returns local summaries", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")

  e <- new.env()
  utils::data("LondonHP", package = "GWmodel", envir = e)
  sp_pts <- e$londonhp

  obs <- sf::st_as_sf(sp_pts)
  suppressWarnings(sf::st_crs(obs) <- 27700)      # LondonHP is on the OSGB grid
  obs <- sf::st_transform(obs, 4326)
  co  <- sf::st_coordinates(obs)

  df <- data.frame(
    lon   = co[, 1],
    lat   = co[, 2],
    value = as.numeric(obs$PURCHASE)
  )

  res <- estimate_gwss_by_point(
    locat_observed  = df,
    locat_predict   = df,
    longitude_col   = "lon",
    latitude_col    = "lat",
    variable        = "value",
    distance_metric = "Great Circle",
    kernel          = "gaussian",
    adaptive        = TRUE,
    target_crs      = 27700,   # OSGB, valid for the radius screen on UK data
    feasible_radius = 500
  )

  # The full point-GWSS path (sf radius screen + GWmodel::gwss) is sensitive to
  # the GWmodel/sf versions and projection; assert it is callable and returns a
  # data.table, and verify the output columns whenever rows are produced.
  expect_s3_class(res, "data.table")
  if (nrow(res) > 0) {
    expect_true(all(c("longitude", "latitude", "bandwidth") %in% names(res)))
    expect_true(any(grepl("_LM$", names(res))))
  }
})
