test_that("avgCIWidthPath returns one positive width per dimension", {
  ci = avgCIWidthPath(ExEstPath_FullFours)
  pi = avgCIWidthPath(ExEstPath_FullFours, P_Int = TRUE)

  expect_named(ci, c("X1", "X2"))
  expect_true(all(ci > 0))
  expect_true(all(pi > ci))                      # prediction intervals add the noise
  expect_true(all(avgCIWidthPath(ExEstPath_FullFours, CL = 0.99) > ci))
})

test_that("avgCIWidthPath matches the band drawn by visualizeEstPath", {
  for(P_Int in c(FALSE, TRUE)){
    for(d in c("X1", "X2")){
      band = ggplot2::layer_data(visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "t", y_var = d, P_Int = P_Int), 1)
      expect_equal(unname(avgCIWidthPath(ExEstPath_FullFours, P_Int = P_Int)[d]), mean(band$ymax - band$ymin), tolerance = 1e-8)
    }
  }
})

test_that("avgCIWidthTraj returns one positive width per trajectory", {
  ci = avgCIWidthTraj(ExEstTraj_OneShot)
  expect_named(ci, c("Trajectory1", "Trajectory2"))
  expect_true(all(ci > 0))
  expect_true(all(avgCIWidthTraj(ExEstTraj_OneShot, CL = 0.99) > ci))
})

test_that("avgCIWidthTraj matches the bands drawn by visualizeEstTraj", {
  p = visualizeEstTraj(EstTraj = ExEstTraj_OneShot, EstPath = ExEstPath_FullFours, baseVectorFields = baseVectorFields,
                       t_grid_size = 0.01, CL = 0.95)
  band = ggplot2::layer_data(p, 3)
  plotted = tapply(band$ymax - band$ymin, band$group, mean)
  expect_equal(unname(avgCIWidthTraj(ExEstTraj_OneShot)), unname(as.numeric(plotted)), tolerance = 1e-8)
})

test_that("trajectory intervals have no prediction (sigma^2) term", {
  expect_false("P_Int" %in% names(formals(avgCIWidthTraj)))
  expect_false("P_Int" %in% names(formals(visualizeEstTraj)))

  # Coverage uses the coefficient covariance only: inflating sigma^2 changes nothing
  inflated = ExEstTraj_OneShot
  inflated$Sigma2 = lapply(inflated$Sigma2, function(s) s * 1e6)
  expect_equal(coverageRateTraj(inflated, ExTraj$Traj, ExTraj$TimeSplits),
               coverageRateTraj(ExEstTraj_OneShot, ExTraj$Traj, ExTraj$TimeSplits))
  expect_equal(avgCIWidthTraj(inflated), avgCIWidthTraj(ExEstTraj_OneShot))
})

test_that("the width functions need a covariance and check their arguments", {
  expect_error(avgCIWidthPath(ExEstPath_ByFours), "covariance")
  expect_error(avgCIWidthTraj(ExEstTraj_ByFours), "covariance")
  expect_error(avgCIWidthPath(ExEstPath_FullFours, CL = 1), "CL")
  expect_error(avgCIWidthTraj(ExEstTraj_OneShot, t_grid_size = 0), "t_grid_size")
})
