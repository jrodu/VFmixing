# A small, quick version of the study
small_run = function(seed, progress = FALSE, B = 3, ...){
  runSimulationStudy_ErrorProcess2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                    traj_nsplits = 20, B = B, progress = progress, ...)
}

test_that("the study draws one SNR and returns its vector fields", {
  run = small_run(1, keep_fits = FALSE)

  expect_true(run$SNR >= 0.1 && run$SNR <= 10)
  expect_equal(run$Results$snr, run$SNR)
  # rotation (y, -x)/r and expansion (x, y)/r at (10, -2)
  r = sqrt(104)
  expect_equal(run$VectorFields(3, c(10, -2)), matrix(c(-2, -10, 10, -2)/r, 2))
  # the returned function does not carry the study's data with it
  expect_lt(length(serialize(run$VectorFields, NULL)), 10000)
})

test_that("every setting in a replicate shares the SNR", {
  run = runSimulationStudy_ErrorProcess2D(seed = 2, error_sd = c(low = 0.001, high = 0.1), n_regions = c(low = 2),
                                          t_end = 4, traj_nsplits = 20, variants = "2x BIC", keep_fits = FALSE)
  expect_equal(nrow(run$Results), 2)
  expect_equal(run$Results$snr, rep(run$SNR, 2))
})

test_that("the true path is a noisy solution whose noise shrinks with the SNR", {
  # A fixed SNR through a one-point range
  lo = runSimulationStudy_ErrorProcess2D(seed = 3, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                         traj_nsplits = 20, log10_snr_range = c(-1, -1), keep_fits = FALSE)
  hi = runSimulationStudy_ErrorProcess2D(seed = 3, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                         traj_nsplits = 20, log10_snr_range = c(3, 3), keep_fits = FALSE)
  expect_equal(lo$SNR, 0.1)
  expect_equal(hi$SNR, 1000)

  # Same seed, so same trajectory and start; the deterministic solution is close to the high-SNR path
  d = hi$Data[[1]]
  rk = RungeKutta(startTime = 0, startPos = unlist(d$TrueData[1, -1]), baseVectorFields = hi$VectorFields,
                  TrajList = d$TrueTraj$Traj, TrajTimeSplits = d$TrueTraj$TimeSplits, endTime = 4, t_step = 0.01)
  dist_hi = sqrt(rowSums((as.matrix(d$TrueData[, -1]) - as.matrix(rk[, -1]))^2))
  dist_lo = sqrt(rowSums((as.matrix(lo$Data[[1]]$TrueData[, -1]) - as.matrix(rk[, -1]))^2))
  expect_lt(max(dist_hi), 0.05)
  expect_gt(max(dist_lo), 10 * max(dist_hi))
})

test_that("one replicate records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "SNR", "VectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 1)
  expect_true(is.na(r$error))
  expect_true(all(c("snr", "path_rmse", "traj_avg_sq_dist", "path_coverage_X1", "traj_coverage_Trajectory2",
                    "path_width_X2", "traj_width_Trajectory1", "variant", "time_total", "selection_share", "traj_coverage_within_Trajectory1") %in% names(r)))
  expect_true(r$path_rmse > 0 && r$traj_avg_sq_dist >= 0)
  expect_named(run$Fits[[1]][[1]], c("EstPath", "EstTraj", "Variance"))
})

test_that("a seed reproduces the replicate", {
  a = small_run(5, keep_fits = FALSE)
  b = small_run(5, keep_fits = FALSE)
  c = small_run(6, keep_fits = FALSE)

  expect_identical(a$Data, b$Data)
  expect_identical(a$SNR, b$SNR)
  drop_time = function(d) d[, !grepl("^time_", names(d))]
  expect_equal(drop_time(a$Results), drop_time(b$Results))
  expect_false(identical(a$SNR, c$SNR))
  expect_null(a$Fits)
})

test_that("the study checks its arguments", {
  expect_error(runSimulationStudy_ErrorProcess2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_ErrorProcess2D(1, log10_snr_range = c(1, -1)), "log10_snr_range")
  expect_error(runSimulationStudy_ErrorProcess2D(1, patience = 0), "patience")
})

test_that("progress is printed as the replicate runs", {
  msgs = testthat::capture_messages(small_run(1, keep_fits = FALSE, progress = TRUE))
  expect_true(any(grepl("^[[]seed 1 [0-9:]{8}[]] simulating data set 1/1: ", msgs)))
  expect_true(any(grepl("fit 1/1: ", msgs)))
  expect_true(any(grepl("done in [0-9.]+ min", msgs)))

  expect_length(testthat::capture_messages(small_run(1, keep_fits = FALSE)), 0)
})
