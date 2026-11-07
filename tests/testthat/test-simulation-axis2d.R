# A small, quick version of the study
small_run = function(seed, methods = c("ByFours-ByFours", "FullRegression-OneShot", "FullFours-OneShot", "Projection-OneShot"), ...){
  runSimulationStudy_Axis2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                            traj_nsplits = 20, methods = methods, ...)
}

test_that("the study returns its vector fields", {
  run = small_run(1, methods = "FullFours-OneShot", keep_fits = FALSE)
  expect_equal(run$VectorFields(3, c(10, -2)), matrix(c(0, 1, 1, 0), 2))
  # the returned function does not carry the study's data with it
  expect_lt(length(serialize(run$VectorFields, NULL)), 10000)
})

test_that("one replicate runs every method pair and records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "VectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 4)
  expect_true(all(is.na(r$error)))
  expect_true(all(c("path_rmse", "traj_avg_sq_dist", "path_coverage_X1", "traj_coverage_Trajectory2",
                    "path_width_X2", "traj_width_Trajectory1", "time_path", "time_traj", "time_total") %in% names(r)))
  expect_true(all(r$path_rmse > 0 & r$traj_avg_sq_dist >= 0))

  # ByFours has no covariance, so no coverage or width
  bf = r[r$method == "ByFours-ByFours", ]
  expect_true(is.na(bf$path_coverage_X1) && is.na(bf$traj_width_Trajectory1))
  os = r[r$method != "ByFours-ByFours", ]
  expect_true(all(os$path_coverage_X1 > 0.8 & os$traj_width_Trajectory1 > 0))
})

test_that("a seed reproduces the data, independently of the methods", {
  a = small_run(5, methods = "FullFours-OneShot", keep_fits = FALSE)
  b = small_run(5, methods = c("FullFours-OneShot", "ByFours-ByFours"), keep_fits = FALSE)
  c = small_run(6, methods = "FullFours-OneShot", keep_fits = FALSE)

  expect_identical(a$Data, b$Data)
  drop_time = function(d) d[, !grepl("^time_", names(d))]
  expect_equal(drop_time(a$Results), drop_time(b$Results[b$Results$method == "FullFours-OneShot", ]), ignore_attr = TRUE)
  expect_false(identical(a$Data, c$Data))
  expect_null(a$Fits)
})

test_that("the study checks its arguments", {
  expect_error(small_run(1, methods = "ByFours-OneShot"), "covariance")
  expect_error(small_run(1, methods = "Spline-OneShot"), "Unknown method")
  expect_error(runSimulationStudy_Axis2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_Axis2D(1, patience = 0), "patience")
})

