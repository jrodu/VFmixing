#' Axis 2-D Simulation Study
#'
#' @description
#' One replicate of the simulation study with two constant vector fields: the
#' first points along the second axis \eqn{(0, 1)} and the second along the
#' first axis \eqn{(1, 0)}, everywhere and at all times. For every combination
#' of positional error level and true trajectory complexity it simulates a data
#' set and fits each path and trajectory method pair, recording accuracy,
#' coverage, interval width and run time. Run it with seeds 1, 2, ..., N to
#' reproduce the full study.
#'
#' This study is self-contained: everything it needs is defined inside this
#' function, apart from the package's estimation and evaluation functions.
#'
#' @details
#' **Data.** For each complexity, true trajectories with that many equally
#' spaced regions on \eqn{[0, t_{end}]} are drawn with
#' [generateFullTrajectory()] (`sigma`, `k`). The true path is the
#' [RungeKutta()] solution through the trajectory-weighted fields from a start
#' position drawn uniformly on `start_range` in each dimension, with time step
#' `t_step`. The observed data add independent normal noise with the error
#' level's standard deviation to each position. All data sets are simulated
#' (in a fixed order) before any fitting, so the data for a seed do not depend
#' on which methods are run.
#'
#' **Methods.** Each entry of `methods` is `"<path method>-<trajectory method>"`.
#' Path methods: `ByFours`, `FullRegression`, `FullFours`, `Projection`.
#' Trajectory methods: `ByFours`, `OneShot`. Selection uses BIC with C1
#' smoothness, and `patience` is passed to every search.
#'
#' **Metrics.** For each setting and method pair:
#' * `path_rmse`: root mean squared distance between the estimated and true
#'   positions at the data times.
#' * `traj_avg_sq_dist`: [cubicAvgSquareDist()] between the estimated and true
#'   trajectories.
#' * `path_coverage_*`: [coverageRatePath()] of the data by the prediction
#'   intervals, per dimension.
#' * `traj_coverage_*`: [coverageRateTraj()] of the true trajectories by the
#'   confidence intervals, per trajectory.
#' * `path_width_*`: [avgCIWidthPath()] with prediction intervals, per dimension.
#' * `traj_width_*`: [avgCIWidthTraj()], per trajectory.
#' * `time_path`, `time_traj`, `time_total`: elapsed seconds.
#'
#' Coverage and width need a covariance structure, so they are `NA` for
#' ByFours paths and trajectories. If a method fails, its metrics are `NA` and
#' the error message is recorded.
#'
#' @param seed an integer. The random seed for this replicate.
#' @param error_sd a named numeric vector. The positional error standard deviations, named by level.
#' @param n_regions a named numeric vector. The number of regions in the true trajectories, named by complexity level.
#' @param methods a character vector of method pairs. See Details.
#' @param patience NULL or a positive integer. Passed to every path and trajectory search.
#' @param t_end a positive number. The end time of the simulation (it starts at 0).
#' @param t_step a positive number. The time step of the data and of the Runge-Kutta solution.
#' @param path_nsplits a positive integer, or NULL for a quarter of the number of data points.
#' The maximum number of splits for the path searches.
#' @param traj_nsplits a positive integer. The maximum number of splits for the trajectory searches.
#' @param CL a number between 0 and 1. The level of the intervals.
#' @param sigma,k positive numbers. Passed to [generateFullTrajectory()].
#' @param start_range a numeric vector of length 2. The range of the uniform start position.
#' @param keep_fits a boolean. If TRUE, the fitted paths and trajectories are returned as well.
#' @param quiet a boolean. If TRUE, the progress output of the estimation functions is suppressed.
#' @param progress a boolean. If TRUE, a time-stamped message is printed as each data set is
#' simulated and each fit finishes, and the total time at the end.
#'
#' @returns a list with:
#' * `Seed`: the seed.
#' * `Settings`: every setting used.
#' * `VectorFields`: the base vector fields function used, e.g. for [visualizeEstTraj()].
#' * `Results`: a data frame with one row per setting and method pair, holding the metrics.
#' * `Data`: for each setting, the true trajectory (`TrueTraj`), true path (`TrueData`) and observed data (`SimData`).
#' * `Fits`: if `keep_fits`, the estimated path and trajectory for each setting and method pair.
#' @export
#'
#' @examples
#' \donttest{
#' # One replicate with a single, quick setting
#' run = runSimulationStudy_Axis2D(seed = 1, error_sd = c(med = 0.01),
#'                                 n_regions = c(low = 2), t_end = 5,
#'                                 methods = "FullFours-OneShot")
#' run$Results
#'
#' # The full study: one replicate per seed
#' # runs = lapply(1:1000, function(s) runSimulationStudy_Axis2D(seed = s, keep_fits = FALSE))
#' # results = do.call(rbind, lapply(runs, `[[`, "Results"))
#' }
#'
runSimulationStudy_Axis2D = function(seed,
                                     error_sd = c(low = 0.001, med = 0.01, high = 0.1),
                                     n_regions = c(low = 2, med = 10, high = 50),
                                     methods = c("ByFours-ByFours", "FullRegression-OneShot", "FullFours-OneShot", "Projection-OneShot"),
                                     patience = 10, t_end = 20, t_step = 0.01,
                                     path_nsplits = NULL, traj_nsplits = 100, CL = 0.95,
                                     sigma = 0.5, k = 0.5, start_range = c(-100, 100),
                                     keep_fits = TRUE, quiet = TRUE, progress = TRUE){

  # ---- The vector fields: (0, 1) and (1, 0) everywhere ----------------------

  vectorFields = function(t, curPos){
    matrix(c(0, 1,
             1, 0), nrow = 2)
  }
  # No source references and a clean environment, so saving the returned
  # function doesn't also save the source file or everything else inside
  # this study.
  vectorFields = utils::removeSource(vectorFields)
  environment(vectorFields) = baseenv()

  n_dim = 2
  n_models = 2

  # ---- Check the arguments --------------------------------------------------

  if(is.null(names(error_sd)) || any(names(error_sd) == "") || any(error_sd < 0)){
    stop("error_sd must be a named vector of non-negative standard deviations.", call. = FALSE)
  }

  if(is.null(names(n_regions)) || any(names(n_regions) == "") || any(n_regions < 1) || any(n_regions != round(n_regions))){
    stop("n_regions must be a named vector of positive whole numbers.", call. = FALSE)
  }

  method_pieces = strsplit(methods, "-", fixed = TRUE)
  known = vapply(method_pieces, function(p){
    length(p) == 2 && p[1] %in% c("ByFours", "FullRegression", "FullFours", "Projection") && p[2] %in% c("ByFours", "OneShot")
  }, logical(1))

  if(!all(known)){
    stop("Unknown method pair(s): ", paste(methods[!known], collapse = ", "),
         ". Use \"<ByFours|FullRegression|FullFours|Projection>-<ByFours|OneShot>\".", call. = FALSE)
  }

  if(any(vapply(method_pieces, function(p){p[1] == "ByFours" && p[2] == "OneShot"}, logical(1)))){
    stop("OneShot trajectories need a path with a covariance structure, so they can't follow ByFours.", call. = FALSE)
  }

  if(!is.null(patience) && (length(patience) != 1 || !is.numeric(patience) || patience < 1 || patience != round(patience))){
    stop("patience must be NULL or a positive whole number.", call. = FALSE)
  }

  # ---- Settings --------------------------------------------------------------

  settings = expand.grid(error_level = names(error_sd), complexity = names(n_regions), stringsAsFactors = FALSE)
  settings$error_sd = unname(error_sd[settings$error_level])
  settings$n_true_regions = unname(n_regions[settings$complexity])
  settings$name = stringr::str_c("error_", settings$error_level, ".complexity_", settings$complexity)

  # ---- Simulate one data set -------------------------------------------------

  simulateData = function(n_true_regions, sd){

    traj_splits = seq(0, t_end, length.out = n_true_regions + 1)
    true_traj = generateFullTrajectory(TimeSplits = traj_splits, numModels = n_models, sigma = sigma, k = k)

    start_pos = stats::runif(n_dim, min = start_range[1], max = start_range[2])

    true_data = RungeKutta(startTime = 0, startPos = start_pos, baseVectorFields = vectorFields,
                           TrajList = true_traj, TrajTimeSplits = traj_splits, endTime = t_end, t_step = t_step)

    sim_data = true_data
    for(d in seq_len(n_dim)){
      sim_data[, d + 1] = sim_data[, d + 1] + stats::rnorm(nrow(sim_data), mean = 0, sd = sd)
    }

    list(TrueTraj = list(Traj = true_traj, TimeSplits = traj_splits), TrueData = true_data, SimData = sim_data)

  }

  # ---- Fit one method pair and compute its metrics ---------------------------

  runQuietly = function(expr){
    if(quiet){
      res = NULL
      utils::capture.output(res <- expr)
      res
    } else{
      expr
    }
  }

  orNA = function(expr, default){
    tryCatch(expr, error = function(e){default})
  }

  fitMethodPair = function(method, sim){

    pieces = strsplit(method, "-", fixed = TRUE)[[1]]
    path_method = pieces[1]
    traj_method = pieces[2]

    SimData = sim$SimData
    TrueData = sim$TrueData
    TrueTraj = sim$TrueTraj

    nsplits_path = if(is.null(path_nsplits)) floor(nrow(SimData)/4) else path_nsplits

    EstPath = NULL
    EstTraj = NULL
    error = NA_character_
    time_path = NA_real_
    time_traj = NA_real_

    tryCatch({

      t0 = Sys.time()
      EstPath = runQuietly(switch(path_method,
        ByFours = getEstimatedPaths_ByFours(SimData, nsplits = nsplits_path, Model_Selection_Type = "BIC", patience = patience),
        FullRegression = getEstimatedPaths_FullRegression(SimData, nsplits = nsplits_path, V_Smooth = F, Model_Selection_Type = "BIC", patience = patience),
        FullFours = getEstimatedPaths_FullFours(SimData, nsplits = nsplits_path, V_Smooth = F, Model_Selection_Type = "BIC", patience = patience),
        Projection = getEstimatedPaths_Projection(SimData, nsplits = nsplits_path, V_Smooth = F, Model_Selection_Type = "BIC", patience = patience)))
      time_path = as.numeric(difftime(Sys.time(), t0, units = "secs"))

      t0 = Sys.time()
      EstTraj = runQuietly(switch(traj_method,
        ByFours = getEstimatedTrajectories_ByFours(EstPath = EstPath, nsplits = traj_nsplits, baseVectorFields = vectorFields,
                                                   data_int = t_step, patience = patience),
        OneShot = getEstimatedTrajectories_OneShot(EstPath = EstPath, nsplits = traj_nsplits, baseVectorFields = vectorFields,
                                                   V_Smooth = F, Model_Selection_Type = "BIC", patience = patience)))
      time_traj = as.numeric(difftime(Sys.time(), t0, units = "secs"))

    }, error = function(e){
      error <<- conditionMessage(e)
    })

    dim_names = stringr::str_c("X", seq_len(n_dim))
    traj_names = stringr::str_c("Trajectory", seq_len(n_models))
    na_dims = stats::setNames(rep(NA_real_, n_dim), dim_names)
    na_trajs = stats::setNames(rep(NA_real_, n_models), traj_names)

    path_rmse = NA_real_
    path_cov = path_width = na_dims
    traj_dist = NA_real_
    traj_cov = traj_width = na_trajs

    if(!is.null(EstPath)){
      est_pos = getEstPathPosition(TrueData$t, EstPath)
      path_rmse = sqrt(mean(rowSums((as.matrix(TrueData[, -1]) - est_pos)^2)))

      if(!is.null(EstPath$Cov)){
        path_cov = orNA(stats::setNames(coverageRatePath(EstPath, SimData, CL = CL), dim_names), na_dims)
        path_width = orNA(stats::setNames(avgCIWidthPath(EstPath, CL = CL, P_Int = TRUE, t_grid_size = t_step), dim_names), na_dims)
      }
    }

    if(!is.null(EstTraj)){
      traj_dist = orNA(cubicAvgSquareDist(EstTraj$Traj, TrueTraj$Traj, EstTraj$TimeSplits, TrueTraj$TimeSplits), NA_real_)

      if(!is.null(EstTraj$Cov)){
        traj_cov = orNA(stats::setNames(coverageRateTraj(EstTraj, TrueTraj$Traj, TrueTraj$TimeSplits, CL = CL), traj_names), na_trajs)
        traj_width = orNA(stats::setNames(avgCIWidthTraj(EstTraj, CL = CL, t_grid_size = t_step), traj_names), na_trajs)
      }
    }

    metrics = data.frame(
      method = method, path_method = path_method, traj_method = traj_method,
      path_regions = if(is.null(EstPath)) NA_integer_ else length(EstPath$Path),
      traj_regions = if(is.null(EstTraj)) NA_integer_ else length(EstTraj$Traj),
      path_rmse = path_rmse, traj_avg_sq_dist = traj_dist,
      t(stats::setNames(path_cov, stringr::str_c("path_coverage_", dim_names))),
      t(stats::setNames(traj_cov, stringr::str_c("traj_coverage_", traj_names))),
      t(stats::setNames(path_width, stringr::str_c("path_width_", dim_names))),
      t(stats::setNames(traj_width, stringr::str_c("traj_width_", traj_names))),
      time_path = time_path, time_traj = time_traj, time_total = time_path + time_traj,
      error = error, check.names = FALSE, stringsAsFactors = FALSE)

    list(metrics = metrics, EstPath = EstPath, EstTraj = EstTraj)

  }

  # ---- Run the replicate -----------------------------------------------------

  # Progress messages (to stderr, so `quiet` doesn't hide them)
  say = function(...){
    if(progress){
      message("[seed ", seed, " ", format(Sys.time(), "%H:%M:%S"), "] ", ...)
    }
  }
  run_start = Sys.time()

  # Simulate every data set first so the data for a seed never depend on the methods.
  set.seed(seed)
  sim_data = lapply(seq_len(nrow(settings)), function(i){
    say("simulating data set ", i, "/", nrow(settings), ": ", settings$name[i])
    simulateData(n_true_regions = settings$n_true_regions[i], sd = settings$error_sd[i])
  })
  names(sim_data) = settings$name

  results = list()
  fits = list()

  for(i in seq_len(nrow(settings))){

    fits[[settings$name[i]]] = list()

    for(method in methods){

      fit = fitMethodPair(method, sim_data[[i]])
      say("fit ", length(results) + 1, "/", nrow(settings)*length(methods), ": ", settings$name[i], ", ", method,
          sprintf(" (%.1f s)", fit$metrics$time_total),
          if(!is.na(fit$metrics$error)) stringr::str_c(" FAILED: ", fit$metrics$error))

      results[[length(results) + 1]] = data.frame(seed = seed, settings[i, c("error_level", "error_sd", "complexity", "n_true_regions")],
                                                  fit$metrics, row.names = NULL, check.names = FALSE)

      if(keep_fits){
        fits[[settings$name[i]]][[method]] = list(EstPath = fit$EstPath, EstTraj = fit$EstTraj)
      }
    }
  }

  say("done in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), run_start, units = "mins"))), " min")

  out = list(Seed = seed,
             Settings = list(study = "Axis2D", error_sd = error_sd, n_regions = n_regions, methods = methods,
                             patience = patience, t_end = t_end, t_step = t_step,
                             path_nsplits = path_nsplits, traj_nsplits = traj_nsplits, CL = CL,
                             sigma = sigma, k = k, start_range = start_range),
             VectorFields = vectorFields,
             Results = do.call(rbind, results),
             Data = sim_data)

  if(keep_fits){
    out$Fits = fits
  }

  out

}
