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

  # The interval is the first layer, so it sits under the data
  band = function(p){ d = ggplot2::layer_data(p, 1); mean(d$ymax - d$ymin) }
  expect_gt(band(ci_plot), 0)
  expect_gt(band(pi_plot), band(ci_plot))   # the prediction band adds the observation noise

  # The plotted 95% prediction band contains the same share of data as coverageRatePath()
  d = ggplot2::layer_data(pi_plot, 1)
  inside = mean(ExData$X1 >= stats::approx(d$x, d$ymin, ExData$t)$y & ExData$X1 <= stats::approx(d$x, d$ymax, ExData$t)$y)
  expect_equal(inside, coverageRatePath(ExEstPath_FullFours, ExData)[1], tolerance = 0.01)
})

test_that("visualizeEstPath draws a swept joint region for position plots", {
  p = visualizeEstPath(ExData, ExEstPath_FullFours, t_grid_size = 0.1, x_var = "X1", y_var = "X2")
  region = ggplot2::layer_data(p, 1)
  # one joint region per plotted time
  expect_equal(length(unique(region$group)), length(seq(0, 5, by = 0.1)))
  expect_length(p$layers, 3)
  expect_equal(c(p$labels$x, p$labels$y), c("X1", "X2"))
})

test_that("visualizeEstPath plots paths without a covariance and checks its inputs", {
  p = visualizeEstPath(ExData, ExEstPath_ByFours, x_var = "X1", y_var = "X2")
  expect_length(p$layers, 2)   # data and path only
  expect_error(visualizeEstPath(ExData, ExEstPath_FullFours, x_var = "X1", y_var = "X1"), "two different")
  expect_error(visualizeEstPath(ExData, ExEstPath_FullFours, CL = 1.5), "CL")
})

test_that("getTrajIntervals gives the gridded estimates and the intervals the other functions use", {
  out = getTrajIntervals(ExEstTraj_OneShot, t_grid_size = 0.01, labels = c("A", "B"))
  expect_named(out, c("A", "B"))
  expect_named(out$A, c("t", "estimate", "lower", "upper"))
  expect_true(all(out$A$lower <= out$A$estimate & out$A$estimate <= out$A$upper))

  # The same widths as avgCIWidthTraj(), and the same coverage as coverageRateTraj()
  expect_equal(unname(c(mean(out$A$upper - out$A$lower), mean(out$B$upper - out$B$lower))),
               unname(avgCIWidthTraj(ExEstTraj_OneShot, t_grid_size = 0.01)), tolerance = 1e-8)
  truth = VFmixing:::evaluateTrajAt(out$A$t, ExTraj)
  inside = c(mean(truth[, 1] > out$A$lower & truth[, 1] < out$A$upper), mean(truth[, 2] > out$B$lower & truth[, 2] < out$B$upper))
  expect_equal(inside, unname(coverageRateTraj(ExEstTraj_OneShot, ExTraj$Traj, ExTraj$TimeSplits)), tolerance = 1e-8)

  # Default names, NA bounds without a covariance, and label checks
  expect_named(getTrajIntervals(ExEstTraj_OneShot, t_grid_size = 0.5), c("Trajectory1", "Trajectory2"))
  no_cov = getTrajIntervals(ExEstTraj_ByFours, t_grid_size = 0.5)
  expect_true(all(is.na(no_cov[[1]]$lower)))
  expect_error(getTrajIntervals(ExEstTraj_OneShot, labels = "only one"), "labels")
})

test_that("visualizeEstTraj works without the path and fields, and uses the labels", {
  p = visualizeEstTraj(ExEstTraj_OneShot, t_grid_size = 0.05, labels = c("Geostrophic", "Ekman"))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$layers[[1]]$data$Model), c("Geostrophic", "Ekman"))
  # no collinearity shading without the path and fields: path, ribbon
  expect_length(p$layers, 2)

  # with the path, fields and truth: path, shading, ribbon, truth
  full = visualizeEstTraj(ExEstTraj_OneShot, ExTraj$Traj, ExTraj$TimeSplits, ExEstPath_FullFours, baseVectorFields, t_grid_size = 0.05)
  expect_length(full$layers, 4)
  expect_error(visualizeEstTraj(ExEstTraj_OneShot, TrueTraj = ExTraj$Traj), "TrueTrajTimeSplits")
})
