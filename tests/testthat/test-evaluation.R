test_that("coverageRateTraj returns one rate per model", {
  cov = coverageRateTraj(ExEstTraj_OneShot, ExTraj$Traj, ExTraj$TimeSplits, CL = 0.95)
  expect_length(cov, nrow(ExEstTraj_OneShot$Traj[[1]]))
  expect_true(all(cov >= 0 & cov <= 1))
})

test_that("coverageRatePath returns one rate per dimension", {
  cov = coverageRatePath(ExEstPath_FullFours, ExData, CL = 0.95)
  expect_length(cov, ncol(ExData) - 1)
  expect_true(all(cov > 0.8 & cov <= 1))
})

test_that("cubicAvgSquareDist is zero for identical splines", {
  expect_equal(cubicAvgSquareDist(ExTraj$Traj, ExTraj$Traj, ExTraj$TimeSplits, ExTraj$TimeSplits), 0)
})
