# A small, quick version of the study
small_run = function(seed, progress = FALSE, B = 3, ...){
  runSimulationStudy_Collinearity2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                    traj_nsplits = 20, B = B, progress = progress, ...)
}

test_that("the study draws one angle and returns its vector fields", {
  run = small_run(1, keep_fits = FALSE)

  expect_true(run$Theta > 0 && run$Theta < pi/2)
  expect_equal(run$Results$theta, run$Theta)
  expect_equal(run$VectorFields(3, c(10, -2)), matrix(c(1, 0, cos(run$Theta), sin(run$Theta)), 2))
  # the returned function does not carry the study's data with it
  expect_lt(length(serialize(run$VectorFields, NULL)), 10000)
})

test_that("every setting in a replicate shares the angle", {
  run = runSimulationStudy_Collinearity2D(seed = 2, error_sd = c(low = 0.001, high = 0.1), n_regions = c(low = 2),
                                          t_end = 4, traj_nsplits = 20, variants = "2x BIC", keep_fits = FALSE)
  expect_equal(nrow(run$Results), 2)
  expect_equal(run$Results$theta, rep(run$Theta, 2))
})

test_that("one replicate records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "Theta", "VectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 1)
  expect_true(is.na(r$error))
  expect_true(all(c("theta", "path_rmse", "traj_avg_sq_dist", "path_coverage_X1", "traj_coverage_Trajectory2",
                    "path_width_X2", "traj_width_Trajectory1", "variant", "time_total", "selection_share", "traj_coverage_within_Trajectory1") %in% names(r)))
  expect_true(r$path_rmse > 0 && r$traj_avg_sq_dist >= 0)
  expect_true(r$path_coverage_X1 > 0.8 && r$traj_width_Trajectory1 > 0)
  expect_named(run$Fits[[1]][[1]], c("EstPath", "EstTraj", "Variance"))
})

test_that("a seed reproduces the replicate", {
  a = small_run(5, keep_fits = FALSE)
  b = small_run(5, keep_fits = FALSE)
  c = small_run(6, keep_fits = FALSE)

  expect_identical(a$Data, b$Data)
  expect_identical(a$Theta, b$Theta)
  drop_time = function(d) d[, !grepl("^time_", names(d))]
  expect_equal(drop_time(a$Results), drop_time(b$Results))
  expect_false(identical(a$Theta, c$Theta))
  expect_null(a$Fits)
})

test_that("the study checks its arguments", {
  expect_error(runSimulationStudy_Collinearity2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_Collinearity2D(1, theta_range = c(0, 2)), "theta_range")
  expect_error(runSimulationStudy_Collinearity2D(1, patience = 0), "patience")
})

test_that("progress is printed as the replicate runs", {
  msgs = testthat::capture_messages(small_run(1, keep_fits = FALSE, progress = TRUE))
  expect_true(any(grepl("^[[]seed 1 [0-9:]{8}[]] simulating data set 1/1: ", msgs)))
  expect_true(any(grepl("fit 1/1: ", msgs)))
  expect_true(any(grepl("done in [0-9.]+ min", msgs)))

  expect_length(testthat::capture_messages(small_run(1, keep_fits = FALSE)), 0)
})
