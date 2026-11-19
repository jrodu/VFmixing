# A small Axis data set: constant fields, 2 true trajectory regions
axisFields = function(t, curPos){matrix(c(0, 1, 1, 0), 2)}

axisData = function(seed = 1, t_end = 6, sd = 0.01){
  set.seed(seed)
  splits = seq(0, t_end, length.out = 3)
  traj = generateFullTrajectory(TimeSplits = splits, numModels = 2, sigma = 0.5, k = 0.5)
  d = RungeKutta(0, c(0, 0), axisFields, traj, splits, t_end, 0.01)
  d[, 2:3] = d[, 2:3] + matrix(stats::rnorm(2*nrow(d), 0, sd), ncol = 2)
  d
}

quietly = function(expr){
  res = NULL
  utils::capture.output(res <- expr)
  res
}

test_that("the undersmoothed choice is the first candidate with twice the regions", {
  expect_equal(undersmoothIndex(c(1, 2, 3, 4, 5, 6), 2), 4)
  expect_equal(undersmoothIndex(c(1, 4, 8, 12), 2), 3)
  # none reaches twice: the last candidate
  expect_equal(undersmoothIndex(c(1, 2, 3), 3), 3)
})

test_that("patience fixes the best candidate, then the search runs on to twice its regions", {
  # values: best at the 3rd candidate (3 regions), then no improvement
  cand = lapply(1:12, function(r){list(Path = vector("list", r))})
  vals = c(5, 4, 3, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3)
  score = function(x){vals[length(x$Path)]}
  regions = function(x){length(x$Path)}

  plain = patienceTracker(2, score, which.min, regions, undersmooth = FALSE)
  stops = vapply(cand, plain$add, logical(1))
  expect_equal(which(stops)[1], 5)          # 2 candidates after the best
  expect_equal(plain$best_index(), 3)

  twice = patienceTracker(2, score, which.min, regions, undersmooth = TRUE)
  stops = vapply(cand, twice$add, logical(1))
  expect_equal(which(stops)[1], 6)          # runs on to 2 x 3 = 6 regions
  expect_equal(twice$best_index(), 3)

  # a candidate with twice the regions already exists when patience runs out: stop at once
  wide = patienceTracker(4, score, which.min, regions, undersmooth = TRUE)
  stops = vapply(cand, wide$add, logical(1))
  expect_equal(which(stops)[1], 7)

  # without patience it never stops early, and the best is the overall best
  none = patienceTracker(NULL, score, which.min, regions, undersmooth = TRUE)
  expect_false(any(vapply(cand, none$add, logical(1))))
  expect_equal(none$best_index(), 3)
})

test_that("every path search doubles the selected regions", {
  d = axisData()
  for(f in list(getEstimatedPaths_FullFours, getEstimatedPaths_FullRegression, getEstimatedPaths_Projection)){
    for(p in list(NULL, 10)){
      a = quietly(f(d, nsplits = 150, V_Smooth = FALSE, patience = p, undersmooth = FALSE))
      b = quietly(f(d, nsplits = 150, V_Smooth = FALSE, patience = p, undersmooth = TRUE))
      expect_gte(length(b$Path), 2*length(a$Path))
      expect_lte(length(b$Path), 2*length(a$Path) + 3)   # ByFours candidates can step past 2K
    }
  }
  a = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10, undersmooth = FALSE))
  b = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10, undersmooth = TRUE))
  expect_gte(length(b$Path), 2*length(a$Path))
})

test_that("both trajectory searches double the selected regions", {
  d = axisData()
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  for(p in list(NULL, 10)){
    a = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = p, undersmooth = FALSE))
    b = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = p, undersmooth = TRUE))
    expect_equal(length(b$Traj), 2*length(a$Traj))
  }
  PB = quietly(getEstimatedPaths_ByFours(d, nsplits = 150, patience = 10))
  a = quietly(getEstimatedTrajectories_ByFours(PB, nsplits = 30, baseVectorFields = axisFields, data_int = 0.01, patience = 10, undersmooth = FALSE))
  b = quietly(getEstimatedTrajectories_ByFours(PB, nsplits = 30, baseVectorFields = axisFields, data_int = 0.01, patience = 10, undersmooth = TRUE))
  expect_equal(length(b$Traj), 2*length(a$Traj))
})

test_that("the full method without the bootstrap matches the separate steps", {
  d = axisData()
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  Tj = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10))

  expect_named(fit, c("EstPath", "EstTraj", "Variance", "Bootstrap"))
  expect_null(fit$Variance)
  expect_equal(fit$EstTraj$Traj, Tj$Traj)
  expect_equal(fit$EstTraj$Cov, Tj$Cov)
})

test_that("the bootstrap adds the selection variance to OneShot's covariance", {
  d = axisData()
  base = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  set.seed(3)
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE, B = 4, quiet = TRUE)

  # same estimates and OneShot covariance, plus the conditional means
  expect_equal(fit$EstTraj$Traj, base$EstTraj$Traj)
  expect_equal(fit$EstTraj$Cov, base$EstTraj$Cov)
  expect_equal(fit$Variance$Within, base$EstTraj$Cov)
  expect_length(fit$EstTraj$SelectionMeans, 4)
  expect_named(fit$EstTraj$SelectionMeans[[1]], c("Traj", "TimeSplits"))

  pw = fit$Variance$Pointwise
  expect_named(pw, c("t", "trajectory", "var_within", "var_selection", "var_total", "selection_share"))
  expect_equal(pw$var_total, pw$var_within + pw$var_selection)
  expect_true(all(pw$var_selection >= 0))
  expect_true(all(pw$selection_share >= 0 & pw$selection_share <= 1))

  # the selection variance is the sample variance (divisor B - 1) of the
  # conditional means at each time, with no projection
  t0 = 2.37
  m = vapply(fit$EstTraj$SelectionMeans, function(x){evaluateTrajAt(t0, x)[1, 1]}, 1)
  expect_equal(trajSelectionVariance(fit$EstTraj, t0)[1, 1], sum((m - mean(m))^2)/(length(m) - 1))

  expect_equal(fit$Bootstrap$B, 4)
  expect_equal(fit$Bootstrap$Successful, 4)
  expect_length(fit$Bootstrap$Nodes, 4)
  expect_named(fit$Bootstrap$Nodes[[1]], c("Path", "Traj"))

  # the residuals are scaled to variance RSS/(n - p)
  p = splineFreeParams(length(fit$EstPath$Path), FALSE)
  expect_equal(fit$Bootstrap$ResidualScale, sqrt(nrow(d)/(nrow(d) - p)))

  # the interval functions add the selection variance, and drop it without the means
  expect_true(all(avgCIWidthTraj(fit$EstTraj) > avgCIWidthTraj(base$EstTraj)))
  within = fit$EstTraj
  within$SelectionMeans = NULL
  expect_equal(avgCIWidthTraj(within), avgCIWidthTraj(base$EstTraj))
  truth = list(Traj = list(matrix(c(0, 0, 0, 0.5, 0, 0, 0, 0.5), 2, 4, byrow = TRUE)), TimeSplits = c(0, 6))
  expect_true(all(coverageRateTraj(fit$EstTraj, truth$Traj, truth$TimeSplits) >=
                    coverageRateTraj(base$EstTraj, truth$Traj, truth$TimeSplits)))
})

test_that("refining splits every selected trajectory region evenly", {
  expect_equal(refineTimeSplits(c(0, 1, 3), 2), c(0, 0.5, 1, 2, 3))
  expect_equal(refineTimeSplits(c(0, 3), 3), c(0, 1, 2, 3))
  expect_equal(refineTimeSplits(c(0, 1, 3), 1), c(0, 1, 3))

  d = axisData()
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  a = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10))
  one = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10, refine = 1))
  b = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10, refine = 2))

  expect_equal(one$Traj, a$Traj)
  expect_equal(b$TimeSplits, refineTimeSplits(a$TimeSplits, 2))
  expect_length(b$Traj, 2*length(a$Traj))
  expect_length(b$Cov, 2*length(a$Traj))
  # Finer regions: more variance
  expect_true(mean(avgCIWidthTraj(b)) > mean(avgCIWidthTraj(a)))
})

test_that("the full method refines the trajectory, in the bootstrap replicates too", {
  d = axisData()
  base = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  set.seed(3)
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, traj_refine = 2,
                                      bootstrap = TRUE, B = 3, quiet = TRUE)

  expect_equal(fit$EstTraj$TimeSplits, refineTimeSplits(base$EstTraj$TimeSplits, 2))
  expect_length(fit$EstTraj$SelectionMeans, 3)
  for(nodes in fit$Bootstrap$Nodes){
    # every replicate's trajectory regions come in refined pairs of equal length
    w = diff(nodes$Traj)
    expect_equal(w[c(TRUE, FALSE)], w[c(FALSE, TRUE)])
  }
})

test_that("the full method checks its arguments", {
  expect_error(getEstimatedTrajectories_Full(axisData(), axisFields, traj_refine = 0), "refine")
  expect_error(getEstimatedTrajectories_Full(axisData(), axisFields, traj_refine = 1.5), "refine")
  d = axisData()
  expect_error(getEstimatedTrajectories_Full(d, axisFields, path_method = "ByFours"), "covariance")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, Model_Selection_Type = "VarPos"), "BIC")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, bootstrap = TRUE, B = 1), "B must")
  expect_error(getEstimatedTrajectories_Full(d, axisFields, patience = 0), "patience")
})

test_that("refining a full fit matches refitting with traj_refine, without new searches", {
  d = axisData()
  set.seed(11)
  fit = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE, B = 3, quiet = TRUE)
  set.seed(11)
  direct = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, traj_refine = 2,
                                         bootstrap = TRUE, B = 3, quiet = TRUE)
  refined = refineEstimatedTrajectories_Full(fit, d, axisFields, refine = 2)

  expect_equal(refined$EstTraj$Traj, direct$EstTraj$Traj)
  expect_equal(refined$EstTraj$Cov, direct$EstTraj$Cov)
  expect_equal(refined$EstTraj$SelectionMeans, direct$EstTraj$SelectionMeans)
  expect_equal(refined$Variance$Pointwise, direct$Variance$Pointwise)
  expect_equal(refined$Bootstrap$Nodes, direct$Bootstrap$Nodes)

  # without a bootstrap, and with refine = 1
  plain = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, quiet = TRUE)
  expect_length(refineEstimatedTrajectories_Full(plain, d, axisFields, refine = 2)$EstTraj$Traj, 2*length(plain$EstTraj$Traj))
  expect_identical(refineEstimatedTrajectories_Full(fit, d, axisFields, refine = 1), fit)
  expect_error(refineEstimatedTrajectories_Full(fit, d, axisFields, refine = 0), "refine")
})

test_that("the bootstrap gives the same result on any number of cores", {
  skip_on_cran()
  # On Windows the workers load the installed package, which a load_all() session may not match
  skip_if(.Platform$OS.type == "windows" && exists(".__DEVTOOLS__", envir = asNamespace("VFmixing")), "workers would use the installed package")
  d = axisData()
  set.seed(21)
  one = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE, B = 4, quiet = TRUE)
  set.seed(21)
  two = getEstimatedTrajectories_Full(d, axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE, B = 4,
                                      cores = 2, quiet = TRUE)

  expect_equal(two$EstTraj$SelectionMeans, one$EstTraj$SelectionMeans)
  expect_equal(two$Bootstrap$Nodes, one$Bootstrap$Nodes)
  expect_equal(two$Variance$Pointwise, one$Variance$Pointwise)
  expect_error(getEstimatedTrajectories_Full(d, axisFields, bootstrap = TRUE, cores = 0), "cores")
})

test_that("the globals a top-level vector field uses are found for the workers", {
  env = globalenv()
  assign("vf_test_scale", 2, envir = env)
  helper = function(x){x*vf_test_scale}
  environment(helper) = env
  assign("vf_test_helper", helper, envir = env)
  f = function(t, curPos){matrix(vf_test_helper(c(0, 1, 1, 0)), 2)}
  environment(f) = env
  on.exit(rm(vf_test_scale, vf_test_helper, envir = env))

  expect_setequal(fieldGlobals(f), c("vf_test_helper", "vf_test_scale"))
  # a field defined inside a function carries its own environment
  expect_length(fieldGlobals(axisFields_local <- local(function(t, curPos){matrix(0, 2, 2)})), 0)
})

test_that("forking falls back to fewer cores, and reruns jobs whose worker died", {
  calls = c()
  testthat::local_mocked_bindings(
    mclapply = function(X, FUN, ..., mc.cores){
      calls <<- c(calls, mc.cores)
      if(mc.cores > 2) stop("unable to fork, possible reason: Cannot allocate memory")
      res = lapply(X, FUN)
      res[2] = list(NULL)   # this worker died
      res
    },
    .package = "parallel")

  w = testthat::capture_warnings(res <- forkWithFallback(1:4, function(x){x^2}, 8))
  expect_equal(sum(grepl("Couldn't start", w)), 2)
  expect_equal(sum(grepl("didn't finish", w)), 1)
  expect_equal(res, list(1, 4, 9, 16))
  expect_equal(calls, c(8, 4, 2))

  # other errors are not retried
  testthat::local_mocked_bindings(mclapply = function(...){stop("something else")}, .package = "parallel")
  expect_error(forkWithFallback(1:4, identity, 4), "something else")
})

test_that("the parallel type is sockets unless forks are requested (never on Windows)", {
  old_opt = options(VFmixing.parallel = NULL)
  on.exit(options(old_opt))

  expect_equal(parallelType(), "socket")
  options(VFmixing.parallel = "fork")
  expect_equal(parallelType(), if(.Platform$OS.type == "windows") "socket" else "fork")
})

test_that("parallel jobs report as they finish, and only when asked", {
  job = parallelJob(c(5, 7), function(x){x + 1}, progress = TRUE, label = "Test job", limit_threads = FALSE)
  expect_message(out <- job(2), "Test job 2/2 done")
  expect_equal(out, 8)
  quiet_job = parallelJob(c(5, 7), function(x){x + 1}, progress = FALSE, label = "Test job", limit_threads = FALSE)
  expect_silent(quiet_job(1))
})

test_that("OneShot can skip the returned fit's covariance", {
  d = axisData()
  P = quietly(getEstimatedPaths_FullFours(d, nsplits = 150, patience = 10))
  a = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10))
  b = quietly(getEstimatedTrajectories_OneShot(P, nsplits = 30, baseVectorFields = axisFields, patience = 10, return_cov = FALSE))
  expect_null(b$Cov)
  expect_equal(b$TimeSplits, a$TimeSplits)
  expect_equal(b$Traj, a$Traj, tolerance = 1e-7)
})

test_that("worker_setup runs in every separate worker before the jobs", {
  skip_on_cran()
  old = options(VFmixing.parallel = "socket")
  on.exit(options(old))
  setup = function(){assign("setup_value", 41, envir = globalenv())}
  res = runInParallel(1:4, function(i){get("setup_value", envir = globalenv()) + 1}, cores = 2, worker_setup = setup)
  expect_equal(unlist(res), rep(42, 4))
})

test_that("a compiled-code failure in the workers says how to fix it", {
  old = options(VFmixing.parallel = "socket")
  on.exit(options(old))
  testthat::local_mocked_bindings(runInParallel = function(X, ...){
    lapply(X, function(i){list(Mean = NULL, Nodes = NULL, Error = "NULL value passed as symbol address")})
  })
  expect_error(getEstimatedTrajectories_Full(axisData(), axisFields, path_nsplits = 150, traj_nsplits = 30, bootstrap = TRUE,
                                             B = 3, cores = 2, quiet = TRUE),
               "worker_setup")
  expect_error(getEstimatedTrajectories_Full(axisData(), axisFields, bootstrap = TRUE, worker_setup = 3), "worker_setup")
})

test_that("the global objects worker_setup uses are copied to the workers too", {
  skip_on_cran()
  old = options(VFmixing.parallel = "socket")
  env = globalenv()
  assign("setup_base_value", 40, envir = env)
  setup = function(){assign("setup_value", setup_base_value + 1, envir = globalenv())}
  environment(setup) = env
  on.exit({options(old); rm(setup_base_value, envir = env)})
  res = runInParallel(1:2, function(i){get("setup_value", envir = globalenv()) + 1}, cores = 2, worker_setup = setup)
  expect_equal(unlist(res), rep(42, 2))
})

test_that("socket workers start single-threaded with this session's libraries, and the session is left as it was", {
  skip_on_cran()
  old = options(VFmixing.parallel = NULL)
  on.exit(options(old))
  before = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "R_LIBS"), unset = NA)

  res = runInParallel(1:2, function(i){c(Sys.getenv("OPENBLAS_NUM_THREADS"), Sys.getenv("OMP_NUM_THREADS"),
                                         Sys.getenv("R_LIBS"))}, cores = 2)
  for(r in res){
    expect_equal(r[1:2], c("1", "1"))
    expect_equal(r[3], paste(.libPaths(), collapse = .Platform$path.sep))
  }
  expect_equal(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "R_LIBS"), unset = NA), before)
})

test_that("setting environment variables temporarily restores them", {
  Sys.setenv(VFMIXING_TEST_A = "x")
  Sys.unsetenv("VFMIXING_TEST_B")
  restore = setEnvTemporarily(c(VFMIXING_TEST_A = "1", VFMIXING_TEST_B = "2"))
  expect_equal(Sys.getenv(c("VFMIXING_TEST_A", "VFMIXING_TEST_B"), names = FALSE), c("1", "2"))
  restore()
  expect_equal(Sys.getenv("VFMIXING_TEST_A"), "x")
  expect_true(is.na(Sys.getenv("VFMIXING_TEST_B", unset = NA)))
  Sys.unsetenv("VFMIXING_TEST_A")
})

test_that("workspace objects are found for field functions defined anywhere", {
  env = globalenv()
  assign("vf_cube_test", matrix(c(0, 1, 1, 0), 2), envir = env)
  assign("vf_lookup_test", function(){vf_cube_test}, envir = env)
  environment(env$vf_lookup_test) = env
  on.exit(rm(vf_cube_test, vf_lookup_test, envir = env))

  # Defined inside a function (a closure), using a workspace object by name
  make_field = function(){
    scale = 1
    function(t, curPos){scale*vf_cube_test}
  }
  environment(make_field) = env
  expect_equal(fieldGlobals(make_field()), "vf_cube_test")

  # A closure with a local helper that uses a workspace helper
  make_field2 = function(){
    helper = function(){vf_lookup_test()}
    function(t, curPos){helper()}
  }
  environment(make_field2) = env
  expect_setequal(fieldGlobals(make_field2()), c("vf_lookup_test", "vf_cube_test"))

  # Package and base functions are not followed
  expect_length(fieldGlobals(stats::median), 0)
})

test_that("objects named in worker_globals are copied to the workers", {
  skip_on_cran()
  old = options(VFmixing.parallel = "socket")
  on.exit(options(old))
  run = function(){
    local_cube = 41   # not in the workspace: found where runInParallel is called
    runInParallel(1:2, function(i){get("local_cube", envir = globalenv()) + 1}, cores = 2,
                  worker_globals = "local_cube", globals_env = environment())
  }
  expect_equal(unlist(run()), c(42, 42))
  expect_error(runInParallel(1:2, identity, cores = 2, worker_globals = "no_such_object_xyz", globals_env = environment()),
               "not found")
  expect_error(getEstimatedTrajectories_Full(axisData(), axisFields, bootstrap = TRUE, worker_globals = 3), "worker_globals")
})
