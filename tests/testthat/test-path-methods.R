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

test_that("every search starts from a single region", {
  bf = quietly(getEstimatedPaths_ByFours(ExData, 10, return_all = TRUE))
  expect_equal(vapply(bf, function(x) length(x$Path), numeric(1))[1:5], 1:5)
  expect_equal(bf[[1]]$TimeSplits, range(ExData$t))

  # Data that one cubic describes well runs through every method
  set.seed(1)
  data = stats::setNames(data.frame(t = seq(0, 1, length.out = 200), X1 = 0, X2 = 0), c("t", "X1", "X2"))
  data$X1 = data$t^3 + stats::rnorm(200, 0, 0.01)
  data$X2 = 2*data$t + stats::rnorm(200, 0, 0.01)
  for(f in list(getEstimatedPaths_ByFours, getEstimatedPaths_FullFours, getEstimatedPaths_FullRegression, getEstimatedPaths_Projection)){
    expect_no_error(quietly(f(data, 10)))
  }
})
