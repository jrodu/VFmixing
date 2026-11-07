# A longer simulated path, so the best model isn't simply the first candidate.
simulated_data = function(){
  set.seed(42)
  traj = generateFullTrajectory(TimeSplits = 0:10, numModels = 2, sigma = 0.1, k = 0.5, startWeights = c(1, 0.3))
  sim = RungeKutta(0, c(1, 1), baseVectorFields, traj, 0:10, endTime = 10, t_step = 0.01)
  sim$X1 = sim$X1 + stats::rnorm(nrow(sim), 0, 0.05)
  sim$X2 = sim$X2 + stats::rnorm(nrow(sim), 0, 0.05)
  sim
}

same_path = function(a, b){
  expect_equal(a$TimeSplits, b$TimeSplits)
  expect_equal(a$Path, b$Path)
}

test_that("patience gives the same path when the best model comes early", {
  sim = simulated_data()

  same_path(quietly(getEstimatedPaths_ByFours(sim, 100, patience = 15)),
            quietly(getEstimatedPaths_ByFours(sim, 100)))
  same_path(quietly(getEstimatedPaths_FullFours(sim, 100, patience = 15)),
            quietly(getEstimatedPaths_FullFours(sim, 100)))
  same_path(quietly(getEstimatedPaths_FullRegression(sim, 40, patience = 15)),
            quietly(getEstimatedPaths_FullRegression(sim, 40)))
  same_path(quietly(getEstimatedPaths_Projection(ExData, 20, patience = 10)),
            quietly(getEstimatedPaths_Projection(ExData, 20)))
})

test_that("patience stops the node search early", {
  # With patience = 2 the search keeps the best candidate plus the two after it.
  all = quietly(getEstimatedPaths_ByFours(ExData, 100, return_all = TRUE))
  best = which.min(VFmixing:::pathSelectionValues(all, "BIC", V_Smooth = FALSE))
  early = quietly(getEstimatedPaths_ByFours(ExData, 100, return_all = TRUE, patience = 2))
  expect_gt(length(all), best + 2)
  expect_length(early, best + 2)
})

test_that("patience works with the covariance-based criteria, including BICVel", {
  for(type in c("VarPosVel", "BICVel")){
    ff = quietly(getEstimatedPaths_FullFours(ExData, 20, Model_Selection_Type = type, patience = 3))
    expect_false(is.null(ff$Cov))
  }
  fr = quietly(getEstimatedPaths_FullRegression(ExData, 20, Model_Selection_Type = "BICVel", patience = 3))
  expect_false(is.null(fr$Cov))
})

test_that("patience is validated by every path method", {
  expect_error(getEstimatedPaths_ByFours(ExData, 10, patience = 0), "patience")
  expect_error(getEstimatedPaths_FullFours(ExData, 10, patience = 1.5), "patience")
  expect_error(getEstimatedPaths_FullRegression(ExData, 10, patience = -1), "patience")
  expect_error(getEstimatedPaths_Projection(ExData, 10, patience = "a"), "patience")
})

test_that("BICVel only considers models before the velocity variance jumps", {
  values = rbind(c(5, 4, 3, 1), c(1, 1.1, 10, 11))   # variance jumps after candidate 2
  expect_equal(VFmixing:::pathBestIndex(values, "BICVel"), 2)
  values[2, 3:4] = c(1.2, 1.3)                         # no jump
  expect_equal(VFmixing:::pathBestIndex(values, "BICVel"), 4)
})
