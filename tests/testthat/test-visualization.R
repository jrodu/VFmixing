test_that("visualizeData plots any pair of columns", {
  expect_s3_class(visualizeData(ExData), "ggplot")
  expect_s3_class(visualizeData(ExData, x_var = "t", y_var = "X1"), "ggplot")
  expect_error(visualizeData(ExData, x_var = "X3"), "must be among")
})

test_that("animatePathTraj rejects paths that aren't two-dimensional", {
  skip_if_not_installed("gganimate")
  path3 = list(matrix(0, nrow = 3, ncol = 4))
  expect_error(animatePathTraj(path3, ExTraj$Traj, c(0, 5), ExTraj$TimeSplits, baseVectorFields, 10), "two-dimensional")
})

test_that("getPropagationList keeps Lon/Lat names for 2-D paths", {
  prop = quietly(getPropagationList(Path = ExEstPath_FullFours, TrajList = list(ExEstTraj_OneShot),
                                    baseVectorFields = baseVectorFields, prop_t_grid_size = 1, n_prop_steps = 2))
  expect_equal(colnames(prop[[1]]), c("Time", "LonM1", "LatM1"))
  expect_equal(nrow(prop[[1]]), 5)
})
