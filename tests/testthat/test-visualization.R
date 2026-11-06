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

test_that("visualizeEstPath draws prediction and confidence bands over time", {
  pi_plot = visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "t", y_var = "X1")
  ci_plot = visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "t", y_var = "X1", P_Int = FALSE)
  expect_s3_class(pi_plot, "ggplot")

  band = function(p){ d = ggplot2::layer_data(p, 2); mean(d$ymax - d$ymin) }
  expect_gt(band(ci_plot), 0)
  expect_gt(band(pi_plot), band(ci_plot))   # the prediction band adds the observation noise

  # The plotted 95% prediction band contains the same share of data as coverageRatePath()
  d = ggplot2::layer_data(pi_plot, 2)
  inside = mean(ExData$X1 >= stats::approx(d$x, d$ymin, ExData$t)$y & ExData$X1 <= stats::approx(d$x, d$ymax, ExData$t)$y)
  expect_equal(inside, coverageRatePath(ExEstPath_FullFours, ExData)[1], tolerance = 0.01)
})

test_that("visualizeEstPath draws ellipses for position plots", {
  p = visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "X1", y_var = "X2", n_ellipses = 7)
  ell = ggplot2::layer_data(p, 2)
  expect_equal(length(unique(ell$group)), 7)
})

test_that("visualizeEstPath plots paths without a covariance and checks its inputs", {
  p = visualizeEstPath(ExData, ExEstPath_ByFours, x_var = "X1", y_var = "X2")
  expect_length(p$layers, 2)   # data and path only
  expect_error(visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "X1", y_var = "X1"), "two different")
  expect_error(visualizeEstPath(ExData, ExEstPath_FullFours, CL = 1.5), "CL")
})
