test_that("getBestSplitsTraj validates type and doesn't need a global n_dim", {
  if(exists("n_dim", envir = globalenv())) rm("n_dim", envir = globalenv())

  expect_error(getBestSplitsTraj(list(ExEstTraj_OneShot), ExEstPath_FullFours, type = "nope"), "must be one of")
  expect_error(getBestSplitsTraj(list(ExEstTraj_OneShot), ExEstPath_FullFours, type = "BICPos"), "no longer available")
  expect_error(getBestSplitsTraj(list(ExEstTraj_ByFours), ExEstPath_ByFours, type = "VarPos"), "covariance")
})

test_that("BIC, BICOS and AIC work for both OneShot and ByFours trajectories", {
  for(type in c("BIC", "BICOS", "BICOpt", "AIC")){
    expect_no_error(quietly(getBestSplitsTraj(list(ExEstTraj_OneShot), ExEstPath_FullFours, type = type, V_Smooth = FALSE)))
    expect_no_error(quietly(getBestSplitsTraj(list(ExEstTraj_ByFours), ExEstPath_ByFours, type = type, V_Smooth = FALSE)))
  }
  expect_no_error(quietly(getBestSplitsTraj(list(ExEstTraj_OneShot), ExEstPath_FullFours, type = "VarPos", V_Smooth = FALSE)))
})

test_that("BIC picks the candidate with the smallest error", {
  candidates = list(ExEstTraj_OneShot, ExEstTraj_OneShot)
  candidates[[2]]$Error = candidates[[2]]$Error * 100
  best = quietly(getBestSplitsTraj(candidates, ExEstPath_FullFours, type = "BIC", V_Smooth = FALSE))
  expect_identical(best, candidates[[1]])
})

test_that("getBestSplitsPath BICVel runs on full-regression candidates", {
  data = stats::setNames(data.frame(ExData), c("t", "X1", "X2"))
  node_sets = list(stats::fivenum(data$t), seq(0, 5, length.out = 7), seq(0, 5, length.out = 9))
  candidates = lapply(node_sets, function(ts) runEstPathStep_FullRegression(data, ts, V_Smooth = FALSE))

  best = quietly(getBestSplitsPath(candidates, show = FALSE, V_Smooth = FALSE, type = "BICVel"))
  expect_true(any(vapply(candidates, function(x) identical(x$TimeSplits, best$TimeSplits), logical(1))))
})

test_that("getBestSplitsPath validates type and needs a covariance for variance criteria", {
  expect_error(getBestSplitsPath(list(ExEstPath_ByFours), type = "nope"), "must be one of")
  expect_error(getBestSplitsPath(list(ExEstPath_ByFours), type = "VarPos"), "covariance")
})
