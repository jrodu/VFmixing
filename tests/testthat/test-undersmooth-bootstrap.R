# A small Axis data set: constant fields, 2 true trajectory regions
axisFields = function(t, curPos){matrix(c(0, 1, 1, 0), 2)}

axisData = function(seed = 1, t_end = 6, sd = 0.01){
  set.seed(seed)
  splits = seq(0, t_end, length.out = 3)
  traj = generateFullTrajectory(TimeSplits = splits, numModels = 2, sigma = 0.5, k = 0.5)
  d = RungeKutta(0, c(0, 0), axisFields, traj, splits, t_end, 0.01)
  d[, 2:3] = d[, 2:3] + matrix(stats::rnorm(2*nrow(d), 0, sd), ncol = 2)
  d
}

quietly = function(expr){
  res = NULL
  utils::capture.output(res <- expr)
  res
}

test_that("the undersmoothed choice is the first candidate with twice the regions", {
  expect_equal(undersmoothIndex(c(1, 2, 3, 4, 5, 6), 2), 4)
  expect_equal(undersmoothIndex(c(1, 4, 8, 12), 2), 3)
  # none reaches twice: the last candidate
  expect_equal(undersmoothIndex(c(1, 2, 3), 3), 3)
})

test_that("patience fixes the best candidate, then the search runs on to twice its regions", {
  # values: best at the 3rd candidate (3 regions), then no improvement
  cand = lapply(1:12, function(r){list(Path = vector("list", r))})
  vals = c(5, 4, 3, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3)
  score = function(x){vals[length(x$Path)]}
  regions = function(x){length(x$Path)}

  plain = patienceTracker(2, score, which.min, regions, undersmooth = FALSE)
  stops = vapply(cand, plain$add, logical(1))
  expect_equal(which(stops)[1], 5)          # 2 candidates after the best
  expect_equal(plain$best_index(), 3)

  twice = patienceTracker(2, score, which.min, regions, undersmooth = TRUE)
  stops = vapply(cand, twice$add, logical(1))
  expect_equal(which(stops)[1], 6)          # runs on to 2 x 3 = 6 regions
  expect_equal(twice$best_index(), 3)

  # a candidate with twice the regions already exists when patience runs out: stop at once
  wide = patienceTracker(4, score, which.min, regions, undersmooth = TRUE)
  stops = vapply(cand, wide$add, logical(1))
  expect_equal(which(stops)[1], 7)

  # without patience it never stops early, and the best is the overall best
  none = patienceTracker(NULL, score, which.min, regions, undersmooth = TRUE)
  expect_false(any(vapply(cand, none$add, logical(1))))
  expect_equal(none$best_index(), 3)
})

test_that("every path search doubles the selected regions", {
  d = axisData()
  for(f in list(getEstimatedPaths_FullFours, getEstimatedPaths_FullRegression, getEstimatedPaths_Projection)){
    for(p in list(NULL, 10)){
      a = quietly(f(d, nsplits = 150, V_Smooth = FALSE, patience = p, undersmooth = FALSE))
      b = quietly(f(d, nsplits = 150, V_Smooth = FALSE, patience = p, undersmooth = TRUE))
      expect_gte(length(b$Path), 2*length(a$Path))
      expect_lte(length(b$Path), 2*length(a$Path) + 3)   # ByFours candidates can step past 2K
    }
  }
  a = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10, undersmooth = FALSE))
  b = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10, undersmooth = TRUE))
  expect_gte(length(b$Path), 2*length(a$Path))
})

test_that("both trajectory searches double the selected regions", {
  d = axisData()
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  for(p in list(NULL, 10)){
    a = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = p, undersmooth = FALSE))
    b = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = p, undersmooth = TRUE))
    expect_equal(length(b$Traj), 2*length(a$Traj))
  }
  PB = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10))
  a = quietly(getEstimatedTrajectories_ByFours(PB, nsplits = 30, baseVectorFields = axisFields, data_int = 0.01, patience = 10, undersmooth = FALSE))
  b = quietly(getEstimatedTrajectories_ByFours(PB, nsplits = 30, baseVectorFields = axisFields, data_int = 0.01, patience = 10, undersmooth = TRUE))
  expect_equal(length(b$Traj), 2*length(a$Traj))
})

test_that("the full method without the bootstrap matches the separate steps", {
  d = axisData()
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  Tj = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10))

  expect_named(fit, c("EstPath", "EstTraj", "Variance", "Bootstrap"))
  expect_null(fit$Variance)
  expect_equal(fit$EstTraj$Traj, Tj$Traj)
  expect_equal(fit$EstTraj$Cov, Tj$Cov)
})

test_that("the bootstrap adds the selection variance to OneShot's covariance", {
  d = axisData()
  base = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  set.seed(3)
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE, B = 4, quiet = TRUE)

  # same estimates, more variance
  expect_equal(fit$EstTraj$Traj, base$EstTraj$Traj)
  expect_equal(fit$Variance$Within, base$EstTraj$Cov)
  expect_equal(fit$EstTraj$Cov, mapply(`+`, fit$Variance$Within, fit$Variance$Selection, SIMPLIFY = FALSE))
  expect_equal(fit$Variance$Total, fit$EstTraj$Cov)
  expect_true(all(vapply(fit$Variance$Selection, function(S){min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)}, 1) > -1e-10))

  pw = fit$Variance$Pointwise
  expect_named(pw, c("t", "trajectory", "var_within", "var_selection", "var_total", "selection_share"))
  expect_equal(pw$var_total, pw$var_within + pw$var_selection)
  expect_true(all(pw$selection_share >= 0 & pw$selection_share <= 1))

  expect_equal(fit$Bootstrap$B, 4)
  expect_equal(fit$Bootstrap$Successful, 4)
  expect_length(fit$Bootstrap$Nodes, 4)
  expect_named(fit$Bootstrap$Nodes[[1]], c("Path", "Traj"))

  # the existing interval functions use the total covariance
  expect_true(all(avgCIWidthTraj(fit$EstTraj) >= avgCIWidthTraj(base$EstTraj)))
})

test_that("the full method checks its arguments", {
  d = axisData()
  expect_error(getEstimatedTrajectories_Full(d, axisFields, path_method = "ByFours"), "covariance")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, Model_Selection_Type = "VarPos"), "BIC")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, bootstrap = TRUE, B = 1), "B must")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, patience = 0), "patience")
})
