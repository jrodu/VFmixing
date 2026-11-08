skip_if_not_installed("akima")

# A small, quick version of the study
small_run = function(seed, progress = FALSE, continuity = c("discontinuous", "continuous", "differentiable"), B = 3, ...){
  runSimulationStudy_Continuity2D(seed = seed, error_sd = c(med = 0.01), n_regions = c(low = 2), t_end = 4,
                                  traj_nsplits = 20, continuity = continuity, B = B, progress = progress, ...)
}

test_that("the study draws two grids and builds a field for each continuity", {
  run = small_run(1, continuity = "discontinuous", keep_fits = FALSE)

  expect_length(run$Grids, 2)
  g = run$Grids[[1]]
  expect_named(g, c("x", "y", "u", "v"))
  expect_equal(nrow(g), 21^2)
  expect_equal(max(sqrt(g$u^2 + g$v^2)), 1)

  run = small_run(1, keep_fits = FALSE)
  fields = run$VectorFields
  expect_named(fields, c("discontinuous", "continuous", "differentiable"))

  g1 = run$Grids[[1]]
  g2 = run$Grids[[2]]
  at = function(g, x, y){i = which(g$x == x & g$y == y); c(g$u[i], g$v[i])}

  # every version matches the grid at a grid point
  for(f in fields){
    expect_equal(f(0, c(3, -2)), matrix(c(at(g1, 3, -2), at(g2, 3, -2)), 2), tolerance = 1e-8)
  }

  # nearest grid point is piecewise constant
  expect_equal(fields$discontinuous(0, c(3.4, -1.7)), fields$discontinuous(0, c(3, -2)))

  # bilinear at the centre of a cell is the mean of its corners
  corners = (at(g1, 3, -2) + at(g1, 4, -2) + at(g1, 3, -1) + at(g1, 4, -1))/4
  expect_equal(fields$continuous(0, c(3.5, -1.5))[, 1], corners)

  # the returned functions carry only their grid
  for(f in fields){
    expect_lt(length(serialize(f, NULL)), 30000)
  }
})

test_that("one replicate runs every continuity and records every metric", {
  run = small_run(1)

  expect_named(run, c("Seed", "Settings", "Grids", "VectorFields", "Results", "Data", "Fits"))
  r = run$Results
  expect_equal(nrow(r), 3)
  expect_equal(r$continuity, c("discontinuous", "continuous", "differentiable"))
  expect_true(all(is.na(r$error)))
  expect_true(all(c("path_attempts", "path_rmse", "traj_avg_sq_dist", "path_coverage_X1", "traj_coverage_Trajectory2",
                    "path_width_X2", "traj_width_Trajectory1", "variant", "time_total", "selection_share", "traj_coverage_within_Trajectory1") %in% names(r)))
  expect_true(all(r$path_attempts >= 1 & r$path_rmse > 0))

  # every true path stays on the grid
  for(d in run$Data){
    expect_true(all(abs(as.matrix(d$TrueData[, -1])) <= 10))
  }
})

test_that("a seed reproduces the replicate", {
  a = small_run(5, continuity = "continuous", keep_fits = FALSE)
  b = small_run(5, continuity = "continuous", keep_fits = FALSE)
  c = small_run(6, continuity = "continuous", keep_fits = FALSE)

  expect_identical(a$Data, b$Data)
  expect_identical(a$Grids, b$Grids)
  drop_time = function(d) d[, !grepl("^time_", names(d))]
  expect_equal(drop_time(a$Results), drop_time(b$Results))
  expect_false(identical(a$Grids, c$Grids))
  expect_null(a$Fits)
})

test_that("the study checks its arguments", {
  expect_error(runSimulationStudy_Continuity2D(1, error_sd = c(0.01)), "named")
  expect_error(runSimulationStudy_Continuity2D(1, continuity = "smooth"), "continuity")
  expect_error(runSimulationStudy_Continuity2D(1, start_range = c(-20, 20)), "start_range")
  expect_error(runSimulationStudy_Continuity2D(1, patience = 0), "patience")
})

test_that("progress is printed as the replicate runs", {
  msgs = testthat::capture_messages(small_run(1, continuity = "continuous", keep_fits = FALSE, progress = TRUE))
  expect_true(any(grepl("^[[]seed 1 [0-9:]{8}[]] simulating data set 1/1: ", msgs)))
  expect_true(any(grepl("fit 1/1: ", msgs)))
  expect_true(any(grepl("done in [0-9.]+ min", msgs)))

  expect_length(testthat::capture_messages(small_run(1, continuity = "continuous", keep_fits = FALSE)), 0)
})
