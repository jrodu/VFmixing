# A small, quick version of the study
small_run = function(seed, progress = FALSE, variants = c("BIC", "2x BIC path only", "2x BIC", "2x BIC + bootstrap"), B = 3, ...){
  runSimulationStudy_Axis2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                            traj_nsplits = 20, variants = variants, B = B, progress = progress, ...)
}

test_that("the study returns its vector fields", {
  run = small_run(1, variants = "2x BIC", keep_fits = FALSE)
  expect_equal(run$VectorFields(3, c(10, -2)), matrix(c(0, 1, 1, 0), 2))
  # the returned function does not carry the study's data with it
  expect_lt(length(serialize(run$VectorFields, NULL)), 10000)
})

test_that("one replicate runs every variant and records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "VectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 4)
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
  expect_error(runSimulationStudy_Axis2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_Axis2D(1, patience = 0), "patience")
})


test_that("progress is printed as the replicate runs", {
  msgs = testthat::capture_messages(small_run(1, variants = "2x BIC", keep_fits = FALSE, progress = TRUE))
  expect_true(any(grepl("^[[]seed 1 [0-9:]{8}[]] simulating data set 1/1: ", msgs)))
  expect_true(any(grepl("fit 1/1: ", msgs)))
  expect_true(any(grepl("done in [0-9.]+ min", msgs)))

  expect_length(testthat::capture_messages(small_run(1, variants = "2x BIC", keep_fits = FALSE)), 0)
})

test_that("the variants differ only where they should", {
  run = small_run(2, keep_fits = FALSE)
  r = run$Results
  rows = function(v) r[r$variant == v, ]

  expect_equal(r$variant, c("BIC", "2x BIC path only", "2x BIC", "2x BIC + bootstrap"))
  # the path is doubled for every variant but BIC
  expect_gte(rows("2x BIC")$path_regions, 2*rows("BIC")$path_regions)
  expect_equal(rows("2x BIC path only")$path_regions, rows("2x BIC")$path_regions)
  # the bootstrap leaves the estimates alone and only adds variance
  expect_equal(rows("2x BIC + bootstrap")$traj_regions, rows("2x BIC")$traj_regions)
  expect_equal(rows("2x BIC + bootstrap")$traj_avg_sq_dist, rows("2x BIC")$traj_avg_sq_dist)
  expect_equal(rows("2x BIC + bootstrap")$traj_width_within_Trajectory1, rows("2x BIC")$traj_width_Trajectory1)
})

test_that("the refined bootstrap variant splits every trajectory region in two", {
  run = small_run(7, variants = c("2x BIC + bootstrap", "2x BIC + refine + bootstrap"), keep_fits = FALSE)
  r = run$Results
  expect_true(all(is.na(r$error)))
  plain = r[r$variant == "2x BIC + bootstrap", ]
  refined = r[r$variant == "2x BIC + refine + bootstrap", ]
  expect_equal(refined$traj_regions, 2*plain$traj_regions)
  expect_false(is.na(refined$selection_share))
  expect_true(refined$traj_width_Trajectory1 >= refined$traj_width_within_Trajectory1 - 1e-12)
})

test_that("traj_refine splits the trajectory regions of every variant", {
  a = small_run(8, variants = "2x BIC", keep_fits = FALSE)
  b = small_run(8, variants = "2x BIC", traj_refine = 2, keep_fits = FALSE)
  expect_equal(b$Results$traj_regions, 2*a$Results$traj_regions)
  expect_equal(b$Settings$traj_refine, 2)
  expect_error(small_run(1, traj_refine = 0), "traj_refine")
})
