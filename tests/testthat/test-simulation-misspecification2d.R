# A small, quick version of the study
small_run = function(seed, progress = FALSE, variants = "2x BIC + bootstrap", B = 3, ...){
  runSimulationStudy_Misspecification2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                        traj_nsplits = 20, variants = variants, B = B, progress = progress, ...)
}

test_that("the study returns its true and misspecified vector fields", {
  run = small_run(1, variants = "2x BIC", keep_fits = FALSE)
  # true: rotation (y, -x)/r and expansion (x, y)/r at (10, -2)
  r = sqrt(104)
  expect_equal(run$VectorFields(3, c(10, -2)), matrix(c(-2, -10, 10, -2)/r, 2))
  # fitted: (3y, -4x)/|(3y, 4x)| and (4x, 3y)/|(4x, 3y)| at (10, -2)
  r_fit = sqrt(36 + 1600)
  expect_equal(run$FitVectorFields(3, c(10, -2)), matrix(c(-6, -40, 40, -6)/r_fit, 2))
  # the returned functions do not carry the study's data with them
  expect_lt(length(serialize(run$VectorFields, NULL)), 10000)
  expect_lt(length(serialize(run$FitVectorFields, NULL)), 10000)
})

test_that("one replicate runs every variant and records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "VectorFields", "FitVectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 1)
  expect_true(all(is.na(r$error)))
  expect_true(all(c("path_rmse", "traj_avg_sq_dist", "path_coverage_X1", "traj_coverage_Trajectory2",
                    "path_width_X2", "traj_width_Trajectory1", "variant", "time_total", "selection_share", "traj_coverage_within_Trajectory1") %in% names(r)))
  expect_true(all(r$path_rmse > 0 & r$traj_avg_sq_dist >= 0))

  boot = r[r$variant == "2x BIC + bootstrap", ]
  expect_equal(nrow(boot), 1)
  expect_true(all(boot$path_coverage_X1 > 0.8 & boot$traj_width_Trajectory1 > 0))
  # the bootstrap only adds variance to OneShot's
  expect_true(all(boot$traj_width_Trajectory1 >= boot$traj_width_within_Trajectory1 - 1e-12))
  expect_true(all(boot$selection_share >= 0 & boot$selection_share <= 1))
  # the other variants have no selection part
  expect_true(all(is.na(r$selection_share[r$variant != "2x BIC + bootstrap"])))
})

test_that("a seed reproduces the data, independently of the variants", {
  a = small_run(5, variants = "2x BIC", keep_fits = FALSE)
  b = small_run(5, variants = c("2x BIC", "BIC"), keep_fits = FALSE)
  c = small_run(6, variants = "2x BIC", keep_fits = FALSE)

  expect_identical(a$Data, b$Data)
  drop_time = function(d) d[, !grepl("^time_", names(d))]
  expect_equal(drop_time(a$Results), drop_time(b$Results[b$Results$variant == "2x BIC", ]), ignore_attr = TRUE)
  expect_false(identical(a$Data, c$Data))
  expect_null(a$Fits)
})

test_that("the study checks its arguments", {
  expect_error(small_run(1, variants = "Spline"), "variants")
  expect_error(small_run(1, B = 1), "B must")
  expect_error(runSimulationStudy_Misspecification2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_Misspecification2D(1, patience = 0), "patience")
})


test_that("progress is printed as the replicate runs", {
  msgs = testthat::capture_messages(small_run(1, variants = "2x BIC", keep_fits = FALSE, progress = TRUE))
  expect_true(any(grepl("^[[]seed 1 [0-9:]{8}[]] simulating data set 1/1: ", msgs)))
  expect_true(any(grepl("fit 1/1: ", msgs)))
  expect_true(any(grepl("done in [0-9.]+ min", msgs)))

  expect_length(testthat::capture_messages(small_run(1, variants = "2x BIC", keep_fits = FALSE)), 0)
})
