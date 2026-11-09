#' Lorenz 3-D Simulation Study
#'
#' @description
#' One replicate of the simulation study with three mutually orthogonal vector
#' fields in 3-D built from the Lorenz system
#' \eqn{L = (10(y - x), x(28 - z) - y, xy - 8z/3)}:
#' * Field 1: the Lorenz field \eqn{L} itself.
#' * Field 2: \eqn{(-L_2, L_1, 0)}, orthogonal to \eqn{L} in the xy-plane.
#' * Field 3: \eqn{(-L_1 L_3, -L_2 L_3, L_1^2 + L_2^2)}, orthogonal to both.
#'
#' Each field is scaled to unit length. Without that, field 3 grows like
#' \eqn{\|L\|^2} and the true paths diverge within a few time steps. The fields
#' are undefined at the Lorenz fixed points, where \eqn{L = 0}: the origin and
#' \eqn{(\pm\sqrt{72}, \pm\sqrt{72}, 27)}. For every combination
#' of positional error level and true trajectory complexity it simulates a data
#' set and fits each variant of the method, recording accuracy,
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
#' on which variants are run.
#'
#' **Variants.** Every variant fits a FullFours path and OneShot trajectories
#' with BIC selection, C1 smoothness and `patience`. Each entry of `variants`
#' is one of:
#' * `"BIC"`: the regions BIC selects, for both the path and the trajectories.
#' * `"2x BIC path only"`: twice BIC's path regions, BIC's trajectory regions.
#' * `"2x BIC"`: twice BIC's regions for both (the default undersmoothing).
#' * `"2x BIC + bootstrap"`: `"2x BIC"`, with the node-selection variance from
#'   `B` bootstrap replicates added to OneShot's covariance, through the law of
#'   total variance ([getEstimatedTrajectories_Full()]).
#'
#' **Metrics.** For each setting and variant:
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
#' * `traj_coverage_within_*`, `traj_width_within_*`: for `"2x BIC + bootstrap"`,
#'   the same with OneShot's covariance alone (without the selection part).
#' * `selection_share`: for `"2x BIC + bootstrap"`, the average share of the
#'   trajectory variance that comes from the node selection.
#' * `time_total`: elapsed seconds for the fit (including the bootstrap).
#'
#' If a fit fails, its metrics are `NA` and the error message is recorded.
#'
#' @param seed an integer. The random seed for this replicate.
#' @param error_sd a named numeric vector. The positional error standard deviations, named by level.
#' @param n_regions a named numeric vector. The number of regions in the true trajectories, named by complexity level.
#' @param variants a character vector of variants. See Details.
#' @param B a whole number of at least 2. The number of bootstrap replicates for `"2x BIC + bootstrap"`.
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
#' * `Results`: a data frame with one row per setting and variant, holding the metrics.
#' * `Data`: for each setting, the true trajectory (`TrueTraj`), true path (`TrueData`) and observed data (`SimData`).
#' * `Fits`: if `keep_fits`, the estimated path, trajectory and (with the bootstrap) variance pieces for each setting and variant.
#' @export
#'
#' @examples
#' \donttest{
#' # One replicate with a single, quick setting
#' run = runSimulationStudy_Lorenz3D(seed = 1, error_sd = c(med = 0.01),
#'                                   n_regions = c(low = 2), t_end = 5, variants = "2x BIC")
#' run$Results
#'
#' # The full study: one replicate per seed
#' # runs = lapply(1:100, function(s){
#' #   runSimulationStudy_Lorenz3D(seed = s, keep_fits = FALSE)
#' # })
#' # results = do.call(rbind, lapply(runs, `[[`, "Results"))
#' }
#'
runSimulationStudy_Lorenz3D = function(seed,
                                       error_sd = c(low = 0.001, med = 0.01, high = 0.1),
                                       n_regions = c(low = 2, med = 10, high = 50),
                                       variants = "2x BIC + bootstrap", B = 20,
                                       patience = 10, t_end = 20, t_step = 0.01,
                                       path_nsplits = NULL, traj_nsplits = 100, CL = 0.95,
                                       sigma = 0.5, k = 0.5, start_range = c(-100, 100),
                                       keep_fits = TRUE, quiet = TRUE, progress = TRUE){

  # ---- The vector fields: Lorenz and two orthogonal companions --------------

  vectorFields = function(t, curPos){

    x = curPos[1]
    y = curPos[2]
    z = curPos[3]

    # Field 1: the standard Lorenz field
    lorenz = c(10 * (y - x),
               x * (28 - z) - y,
               x * y - (8/3) * z)

    # Field 2: orthogonal to field 1 in the xy-plane
    rot = c(-lorenz[2],
            lorenz[1],
            0)

    # Field 3: orthogonal to fields 1 and 2
    ortho = c(-lorenz[1] * lorenz[3],
              -lorenz[2] * lorenz[3],
              lorenz[1]^2 + lorenz[2]^2)

    fields = cbind(lorenz, rot, ortho, deparse.level = 0)

    # Scale each field to unit length
    sweep(fields, 2, sqrt(colSums(fields^2)), "/")

  }
  # No source references and a clean environment, so saving the returned
  # function doesn't also save the source file or everything else inside
  # this study.
  vectorFields = utils::removeSource(vectorFields)
  environment(vectorFields) = baseenv()

  n_dim = 3
  n_models = 3

  # ---- Check the arguments --------------------------------------------------

  if(is.null(names(error_sd)) || any(names(error_sd) == "") || any(error_sd < 0)){
    stop("error_sd must be a named vector of non-negative standard deviations.", call. = FALSE)
  }

  if(is.null(names(n_regions)) || any(names(n_regions) == "") || any(n_regions < 1) || any(n_regions != round(n_regions))){
    stop("n_regions must be a named vector of positive whole numbers.", call. = FALSE)
  }

  if(length(variants) < 1 || anyDuplicated(variants) || !all(variants %in% c("BIC", "2x BIC path only", "2x BIC", "2x BIC + bootstrap"))){
    stop("variants must be distinct values from \"BIC\", \"2x BIC path only\", \"2x BIC\" and \"2x BIC + bootstrap\".", call. = FALSE)
  }

  if("2x BIC + bootstrap" %in% variants && (length(B) != 1 || !is.numeric(B) || B < 2 || B != round(B))){
    stop("B must be a whole number of at least 2.", call. = FALSE)
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

  # ---- Fit one variant and compute its metrics -------------------------------

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

  # Every variant fits a FullFours path and OneShot trajectories with BIC
  # selection and C1 smoothness. They differ in the number of regions and in
  # the trajectory covariance:
  # * "BIC": the regions BIC selects, for the path and the trajectories.
  # * "2x BIC path only": twice BIC's path regions, BIC's trajectory regions.
  # * "2x BIC": twice BIC's regions for both.
  # * "2x BIC + bootstrap": "2x BIC", with the node-selection variance added to
  #   OneShot's covariance (getEstimatedTrajectories_Full()).
  fitVariant = function(variant, sim, fitFields, extra = list()){

    SimData = sim$SimData
    TrueData = sim$TrueData
    TrueTraj = sim$TrueTraj

    nsplits_path = if(is.null(path_nsplits)) floor(nrow(SimData)/4) else path_nsplits

    EstPath = NULL
    EstTraj = NULL
    Variance = NULL
    error = NA_character_
    time_total = NA_real_

    tryCatch({

      t0 = Sys.time()

      if(variant == "2x BIC + bootstrap"){
        fit = runQuietly(getEstimatedTrajectories_Full(SimData, fitFields, path_nsplits = nsplits_path, traj_nsplits = traj_nsplits,
                                                       patience = patience, undersmooth = TRUE, bootstrap = TRUE, B = B, quiet = TRUE))
        EstPath = fit$EstPath
        EstTraj = fit$EstTraj
        Variance = fit$Variance
      } else{
        EstPath = runQuietly(getEstimatedPaths_FullFours(SimData, nsplits = nsplits_path, V_Smooth = F, Model_Selection_Type = "BIC",
                                                         patience = patience, undersmooth = variant != "BIC"))
        EstTraj = runQuietly(getEstimatedTrajectories_OneShot(EstPath = EstPath, nsplits = traj_nsplits, baseVectorFields = fitFields,
                                                              V_Smooth = F, Model_Selection_Type = "BIC", patience = patience,
                                                              undersmooth = variant == "2x BIC"))
      }

      time_total = as.numeric(difftime(Sys.time(), t0, units = "secs"))

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
    traj_cov_within = traj_width_within = na_trajs
    selection_share = NA_real_

    if(!is.null(EstPath)){
      est_pos = getEstPathPosition(TrueData$t, EstPath)
      path_rmse = sqrt(mean(rowSums((as.matrix(TrueData[, -1]) - est_pos)^2)))

      path_cov = orNA(stats::setNames(coverageRatePath(EstPath, SimData, CL = CL), dim_names), na_dims)
      path_width = orNA(stats::setNames(avgCIWidthPath(EstPath, CL = CL, P_Int = TRUE, t_grid_size = t_step), dim_names), na_dims)
    }

    if(!is.null(EstTraj)){
      traj_dist = orNA(cubicAvgSquareDist(EstTraj$Traj, TrueTraj$Traj, EstTraj$TimeSplits, TrueTraj$TimeSplits), NA_real_)

      # With the bootstrap, the intervals add the selection variance from EstTraj$SelectionMeans
      traj_cov = orNA(stats::setNames(coverageRateTraj(EstTraj, TrueTraj$Traj, TrueTraj$TimeSplits, CL = CL), traj_names), na_trajs)
      traj_width = orNA(stats::setNames(avgCIWidthTraj(EstTraj, CL = CL, t_grid_size = t_step), traj_names), na_trajs)

      if(!is.null(Variance)){
        # The same fit with OneShot's (within-node) covariance alone
        EstTraj_within = EstTraj
        EstTraj_within$SelectionMeans = NULL
        traj_cov_within = orNA(stats::setNames(coverageRateTraj(EstTraj_within, TrueTraj$Traj, TrueTraj$TimeSplits, CL = CL), traj_names), na_trajs)
        traj_width_within = orNA(stats::setNames(avgCIWidthTraj(EstTraj_within, CL = CL, t_grid_size = t_step), traj_names), na_trajs)
        selection_share = mean(Variance$Pointwise$selection_share)
      }
    }

    metrics = data.frame(
      variant = variant,
      path_regions = if(is.null(EstPath)) NA_integer_ else length(EstPath$Path),
      traj_regions = if(is.null(EstTraj)) NA_integer_ else length(EstTraj$Traj),
      path_rmse = path_rmse, traj_avg_sq_dist = traj_dist,
      t(stats::setNames(path_cov, stringr::str_c("path_coverage_", dim_names))),
      t(stats::setNames(traj_cov, stringr::str_c("traj_coverage_", traj_names))),
      t(stats::setNames(path_width, stringr::str_c("path_width_", dim_names))),
      t(stats::setNames(traj_width, stringr::str_c("traj_width_", traj_names))),
      t(stats::setNames(traj_cov_within, stringr::str_c("traj_coverage_within_", traj_names))),
      t(stats::setNames(traj_width_within, stringr::str_c("traj_width_within_", traj_names))),
      selection_share = selection_share,
      time_total = time_total,
      error = error, check.names = FALSE, stringsAsFactors = FALSE)

    if(length(extra) > 0){
      metrics = data.frame(extra, metrics, check.names = FALSE, stringsAsFactors = FALSE)
    }

    list(metrics = metrics, EstPath = EstPath, EstTraj = EstTraj, Variance = Variance)

  }

  # ---- Run the replicate -----------------------------------------------------

  # Progress messages (to stderr, so `quiet` doesn't hide them)
  say = function(...){
    if(progress){
      message("[seed ", seed, " ", format(Sys.time(), "%H:%M:%S"), "] ", ...)
    }
  }
  run_start = Sys.time()

  # Simulate every data set first so the data for a seed never depend on the variants.
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

    for(variant in variants){

      fit = fitVariant(variant, sim_data[[i]], vectorFields)
      say("fit ", length(results) + 1, "/", nrow(settings)*length(variants), ": ", settings$name[i], ", ", variant,
          sprintf(" (%.1f s)", fit$metrics$time_total),
          if(!is.na(fit$metrics$error)) stringr::str_c(" FAILED: ", fit$metrics$error))

      results[[length(results) + 1]] = data.frame(seed = seed, settings[i, c("error_level", "error_sd", "complexity", "n_true_regions")],
                                                  fit$metrics, row.names = NULL, check.names = FALSE)

      if(keep_fits){
        fits[[settings$name[i]]][[variant]] = list(EstPath = fit$EstPath, EstTraj = fit$EstTraj, Variance = fit$Variance)
      }
    }
  }

  say("done in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), run_start, units = "mins"))), " min")

  out = list(Seed = seed,
             Settings = list(study = "Lorenz3D", error_sd = error_sd, n_regions = n_regions, variants = variants, B = B,
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
