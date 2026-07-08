test_that("gw_distance_metric_names/presets are consistent", {
  nms <- gw_distance_metric_names()
  expect_type(nms, "character")
  expect_true(all(c("Euclidean", "Great Circle") %in% nms))

  presets <- gw_distance_metric_presets()
  expect_type(presets, "list")
  expect_identical(names(presets), nms)
  # every preset is a list(p, theta, longlat)
  expect_true(all(vapply(presets, function(x)
    all(c("p", "theta", "longlat") %in% names(x)), logical(1))))
})

test_that("resolve_distance_metric returns the right shape and errors sensibly", {
  euc <- resolve_distance_metric("Euclidean")
  expect_equal(euc$p, 2)
  expect_false(euc$longlat)

  gc <- resolve_distance_metric("Great Circle")
  expect_true(gc$longlat)

  # unknown metric: NULL when not stopping, error when stopping
  expect_null(resolve_distance_metric("does-not-exist", stop_on_error = FALSE))
  expect_error(resolve_distance_metric("does-not-exist"))
})
