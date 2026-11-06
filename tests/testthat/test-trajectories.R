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
