test_that("ByFours trajectories report regions in n_splits", {
  traj = quietly(getEstimatedTrajectories_ByFours(EstPath = ExEstPath_FullFours, nsplits = 10,
                                                  baseVectorFields = baseVectorFields, data_int = 0.01))
  expect_equal(traj$n_splits, length(traj$Traj))
  expect_equal(nrow(traj$Traj[[1]]), 2)
})

test_that("OneShot works when path and trajectory smoothness differ", {
  path_c2 = quietly(getEstimatedPaths_FullFours(ExData, 10, V_Smooth = TRUE))
  traj = quietly(getEstimatedTrajectories_OneShot(EstPath = path_c2, nsplits = 2,
                                                  baseVectorFields = baseVectorFields, V_Smooth = FALSE))
  expect_equal(traj$n_splits, length(traj$Traj))
  expect_false(is.null(traj$Cov))
})

test_that("ByFours trajectories support patience", {
  full = quietly(getEstimatedTrajectories_ByFours(EstPath = ExEstPath_FullFours, nsplits = 40,
                                                  baseVectorFields = baseVectorFields, data_int = 0.01))
  early = quietly(getEstimatedTrajectories_ByFours(EstPath = ExEstPath_FullFours, nsplits = 40,
                                                   baseVectorFields = baseVectorFields, data_int = 0.01, patience = 10))
  expect_equal(early$TimeSplits, full$TimeSplits)
  expect_equal(early$Traj, full$Traj)

  # patience smaller than the four initial candidates still works
  expect_no_error(quietly(getEstimatedTrajectories_ByFours(EstPath = ExEstPath_FullFours, nsplits = 40,
                                                           baseVectorFields = baseVectorFields, data_int = 0.01, patience = 1)))
  expect_error(getEstimatedTrajectories_ByFours(EstPath = ExEstPath_FullFours, nsplits = 40,
                                                baseVectorFields = baseVectorFields, patience = 0), "patience")
})

test_that("ByFours trajectories accept fields with zero components", {
  path = quietly(getEstimatedPaths_ByFours(ExData, 20))
  expect_no_error(quietly(getEstimatedTrajectories_ByFours(EstPath = path, nsplits = 10,
                                                           baseVectorFields = function(t, curPos){matrix(c(0, 1, 1, 0), 2)}, data_int = 0.01)))
})
