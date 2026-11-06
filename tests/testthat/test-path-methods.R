# Continuity of an estimated path at its interior nodes.
expect_continuous_at_nodes = function(EstPath){
  nodes = EstPath$TimeSplits[-c(1, length(EstPath$TimeSplits))]
  for(i in seq_along(nodes)){
    left = EstPath$Path[[i]]
    right = EstPath$Path[[i + 1]]
    t0 = nodes[i]
    expect_equal(c(evaluateCubic(t0, left, EstPath$TimeSplits[i])), c(evaluateCubic(t0, right, t0)), tolerance = 1e-8)
    expect_equal(c(evaluateCubicVelocity(t0, left, EstPath$TimeSplits[i])), c(evaluateCubicVelocity(t0, right, t0)), tolerance = 1e-8)
  }
}

test_that("n_splits is the number of regions for every path method", {
  bf = quietly(getEstimatedPaths_ByFours(ExData, 20))
  fr = quietly(getEstimatedPaths_FullRegression(ExData, 20, FALSE, "BIC"))
  ff = quietly(getEstimatedPaths_FullFours(ExData, 20))

  expect_equal(bf$n_splits, length(bf$Path))
  expect_equal(fr$n_splits, length(fr$Path))
  expect_equal(ff$n_splits, length(ff$Path))
})

test_that("FullFours returns a FullRegression-style path on the ByFours nodes", {
  bf = quietly(getEstimatedPaths_ByFours(ExData, 20))
  ff = quietly(getEstimatedPaths_FullFours(ExData, 20))
  fr = quietly(getEstimatedPaths_FullRegression(ExData, 20, FALSE, "BIC"))

  expect_equal(ff$TimeSplits, bf$TimeSplits)
  expect_setequal(names(ff), names(fr))
  expect_false(is.null(ff$Cov))
  expect_continuous_at_nodes(ff)
})

test_that("FullFours supports the covariance-based selection criteria", {
  ff = quietly(getEstimatedPaths_FullFours(ExData, 10, Model_Selection_Type = "VarPosVel"))
  expect_false(is.null(ff$Cov))
  expect_equal(ff$n_splits, length(ff$Path))

  expect_error(getEstimatedPaths_FullFours(ExData, 10, Model_Selection_Type = "nope"), "must be one of")
})

test_that("ByFours only accepts criteria that don't need a covariance", {
  expect_error(getEstimatedPaths_ByFours(ExData, 10, Model_Selection_Type = "VarPos"), "BIC")
  expect_no_error(quietly(getEstimatedPaths_ByFours(ExData, 10, Model_Selection_Type = "AIC")))
})

test_that("path smoothness is detected from the free parameters", {
  c1 = quietly(getEstimatedPaths_FullFours(ExData, 10, V_Smooth = FALSE))
  c2 = quietly(getEstimatedPaths_FullFours(ExData, 10, V_Smooth = TRUE))

  expect_false(VFmixing:::inferPathVSmooth(c1))
  expect_true(VFmixing:::inferPathVSmooth(c2))
  expect_false(VFmixing:::inferPathVSmooth(ExEstPath_ByFours))
  expect_continuous_at_nodes(c2)
})
