#' Continuity 2-D Simulation Study
#'
#' @description
#' One replicate of the simulation study on how the smoothness of the vector
#' fields affects estimation. Each replicate draws two random ocean-like
#' fields on a coarse grid over \eqn{[-10, 10]^2}, each made of two eddies
#' and a meandering Gulf Stream-like current, and turns the grid values into
#' vector fields three ways:
#' * `discontinuous`: the value at the nearest grid point (piecewise constant).
#' * `continuous`: bilinear interpolation ([akima::bilinear()]).
#' * `differentiable`: bicubic interpolation ([akima::bicubic()]).
#'
#' Then, for every combination of positional error level, true trajectory
#' complexity and continuity, it simulates a data set and fits FullFours paths
#' followed by OneShot trajectories with the same fields, recording accuracy,
#' coverage, interval width and run time. Run it with seeds 1, 2, ..., N to
#' reproduce the full study. It needs the akima package.
#'
#' This study is self-contained: everything it needs is defined inside this
#' function, apart from the package's estimation and evaluation functions and
#' akima's interpolation.
#'
#' @details
#' **Grids.** Each grid has spacing 1 on \eqn{[-10, 10]^2}. Its value at a
#' grid point is \eqn{10 e_1 + 10 e_2 + 3 g}, where \eqn{e_1, e_2} are eddies
#' (pure rotations about uniformly placed centres, with random direction,
#' decaying with distance) and \eqn{g} is a current that follows a randomly
#' rotated sine curve with random amplitude, frequency, phase and width. The
#' vectors are then scaled so the longest has length 1. The two grids are
#' drawn once for the whole replicate.
#'
#' **Data.** For each setting, true trajectories with the complexity's number
#' of equally spaced regions on \eqn{[0, t_{end}]} are drawn with
#' [generateFullTrajectory()] (`sigma`, `k`). The true path is the
#' [RungeKutta()] solution through the trajectory-weighted fields of the
#' setting's continuity, from a start position drawn uniformly on
#' `start_range` in each dimension, with time step `t_step`. The fields are
#' only defined on the grid, so if the path leaves \eqn{[-10, 10]^2}, the
#' trajectories and start are drawn again (up to `max_attempts` times); the
#' number of draws is recorded. The observed data add independent normal noise
#' with the error level's standard deviation to each position. All data sets
#' are simulated (in a fixed order) before any fitting.
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
#' Each setting is fitted with its own continuity's fields.
#'
#' **Metrics.** For each setting and variant:
#' * `path_attempts`: the number of draws needed for a path inside the grid.
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
#' @param continuity a character vector, any of `"discontinuous"`, `"continuous"` and `"differentiable"`.
#' @param variants a character vector of variants. See Details.
#' @param B a whole number of at least 2. The number of bootstrap replicates for `"2x BIC + bootstrap"`.
#' @param patience NULL or a positive integer. Passed to the path and trajectory searches.
#' @param t_end a positive number. The end time of the simulation (it starts at 0).
#' @param t_step a positive number. The time step of the data and of the Runge-Kutta solution.
#' @param path_nsplits a positive integer, or NULL for a quarter of the number of data points.
#' The maximum number of splits for the path search.
#' @param traj_nsplits a positive integer. The maximum number of splits for the trajectory search.
#' @param CL a number between 0 and 1. The level of the intervals.
#' @param sigma,k positive numbers. Passed to [generateFullTrajectory()].
#' @param start_range a numeric vector of length 2 within \eqn{[-10, 10]}. The range of the uniform start position.
#' @param max_attempts a positive integer. The most draws allowed for a path that stays inside the grid.
#' @param keep_fits a boolean. If TRUE, the fitted paths and trajectories are returned as well.
#' @param quiet a boolean. If TRUE, the progress output of the estimation functions is suppressed.
#' @param progress a boolean. If TRUE, a time-stamped message is printed as each data set is
#' simulated and each fit finishes, and the total time at the end.
#'
#' @returns a list with:
#' * `Seed`: the seed.
#' * `Settings`: every setting used.
#' * `Grids`: the two random grids, data frames with columns `x`, `y`, `u` and `v`.
#' * `VectorFields`: the base vector fields function of each continuity, e.g. for [visualizeEstTraj()].
#' * `Results`: a data frame with one row per setting and variant, holding the metrics.
#' * `Data`: for each setting, the true trajectory (`TrueTraj`), true path (`TrueData`) and observed data (`SimData`).
#' * `Fits`: if `keep_fits`, the estimated path, trajectory and (with the bootstrap) variance pieces for each setting and variant.
#' @export
#'
#' @examples
#' \donttest{
#' if(requireNamespace("akima", quietly = TRUE)){
#'   # One replicate with a single, quick setting
#'   run = runSimulationStudy_Continuity2D(seed = 1, error_sd = c(med = 0.01),
#'                                         n_regions = c(low = 2),
#'                                         continuity = "continuous", t_end = 5, variants = "2x BIC")
#'   run$Results
#' }
#'
#' # The full study: one replicate per seed
#' # runs = lapply(1:100, function(s){
#' #   runSimulationStudy_Continuity2D(seed = s, keep_fits = FALSE)
#' # })
#' # results = do.call(rbind, lapply(runs, `[[`, "Results"))
#' }
#'
runSimulationStudy_Continuity2D = function(seed,
                                           error_sd = c(low = 0.001),
                                           n_regions = c(med = 10),
                                           continuity = c("discontinuous", "continuous", "differentiable"),
                                           variants = "2x BIC + bootstrap", B = 20,
                                           patience = 10, t_end = 20, t_step = 0.01,
                                           path_nsplits = NULL, traj_nsplits = 100, CL = 0.95,
                                           sigma = 0.5, k = 0.5, start_range = c(-5, 5), max_attempts = 100,
                                           keep_fits = TRUE, quiet = TRUE, progress = TRUE){

  n_dim = 2
  n_models = 2
  grid_coords = seq(-10, 10, by = 1)

  # ---- Check the arguments --------------------------------------------------

  if(!requireNamespace("akima", quietly = TRUE)){
    stop("This study needs the akima package. Install it with install.packages(\"akima\").", call. = FALSE)
  }

  if(is.null(names(error_sd)) || any(names(error_sd) == "") || any(error_sd < 0)){
    stop("error_sd must be a named vector of non-negative standard deviations.", call. = FALSE)
  }

  if(is.null(names(n_regions)) || any(names(n_regions) == "") || any(n_regions < 1) || any(n_regions != round(n_regions))){
    stop("n_regions must be a named vector of positive whole numbers.", call. = FALSE)
  }

  if(length(continuity) < 1 || !all(continuity %in% c("discontinuous", "continuous", "differentiable")) || anyDuplicated(continuity)){
    stop("continuity must be distinct values from \"discontinuous\", \"continuous\" and \"differentiable\".", call. = FALSE)
  }

  if(length(start_range) != 2 || start_range[1] > start_range[2] || start_range[1] < -10 || start_range[2] > 10){
    stop("start_range must be an increasing pair within [-10, 10], the grid.", call. = FALSE)
  }

  if(!is.null(patience) && (length(patience) != 1 || !is.numeric(patience) || patience < 1 || patience != round(patience))){
    stop("patience must be NULL or a positive whole number.", call. = FALSE)
  }

  if(length(variants) < 1 || anyDuplicated(variants) || !all(variants %in% c("BIC", "2x BIC path only", "2x BIC", "2x BIC + bootstrap"))){
    stop("variants must be distinct values from \"BIC\", \"2x BIC path only\", \"2x BIC\" and \"2x BIC + bootstrap\".", call. = FALSE)
  }

  if("2x BIC + bootstrap" %in% variants && (length(B) != 1 || !is.numeric(B) || B < 2 || B != round(B))){
    stop("B must be a whole number of at least 2.", call. = FALSE)
  }

  # ---- The random grid: two eddies and a Gulf Stream-like current -----------

  createRandomGrid = function(){

    grid_df = expand.grid(x = grid_coords, y = grid_coords)

    # Random locations and directions for the eddies
    eddy1_loc = c(stats::runif(1, -8, 8), stats::runif(1, -8, 8))
    eddy2_loc = c(stats::runif(1, -8, 8), stats::runif(1, -8, 8))
    eddy1_dir = sample(c(-1, 1), 1)
    eddy2_dir = sample(c(-1, 1), 1)

    # The current: a sine wave with random height, bends, shift and width,
    # rotated by a random angle
    stream_amplitude = stats::runif(1, 1, 5)
    stream_frequency = stats::runif(1, 0.5, 3)
    stream_phase = stats::runif(1, 0, 2*pi)
    stream_width = stats::runif(1, 1, 5)
    stream_angle = stats::runif(1, 0, 2*pi)
    cos_a = cos(stream_angle)
    sin_a = sin(stream_angle)

    eddyInfluence = function(pt, loc, dir){
      to_eddy = pt - loc
      dist_sq = sum(to_eddy^2) + 1e-6
      dir * c(-to_eddy[2], to_eddy[1])/(dist_sq*0.5 + 1)
    }

    u = numeric(nrow(grid_df))
    v = numeric(nrow(grid_df))

    for(i in seq_len(nrow(grid_df))){

      pt = c(grid_df$x[i], grid_df$y[i])

      # Rotate the point into the stream's frame
      rotated_px = pt[1]*cos_a + pt[2]*sin_a
      rotated_py = -pt[1]*sin_a + pt[2]*cos_a

      # Strength falls off with distance from the stream's centre line
      stream_center_y = stream_amplitude*sin(stream_frequency*rotated_px + stream_phase)
      stream_strength = exp(-(rotated_py - stream_center_y)^2/(2*stream_width^2))

      # Flow along the tangent of the sine curve, rotated back
      slope = stream_amplitude*stream_frequency*cos(stream_frequency*rotated_px + stream_phase)
      tangent_rotated = c(1, slope)/sqrt(1 + slope^2)
      tangent = c(tangent_rotated[1]*cos_a - tangent_rotated[2]*sin_a,
                  tangent_rotated[1]*sin_a + tangent_rotated[2]*cos_a)

      final_vec = 10*eddyInfluence(pt, eddy1_loc, eddy1_dir) +
                  10*eddyInfluence(pt, eddy2_loc, eddy2_dir) +
                  3*tangent*stream_strength

      u[i] = final_vec[1]
      v[i] = final_vec[2]
    }

    # Scale so the longest vector has length 1
    max_mag = max(sqrt(u^2 + v^2))
    grid_df$u = u/max_mag
    grid_df$v = v/max_mag

    grid_df

  }

  # ---- The vector fields of each continuity ----------------------------------

  # Builds the field function from the grid values. Its environment is set to
  # the base environment below, so the returned function holds only the grid,
  # and saving it doesn't also save the rest of the study.
  fieldFactory = function(type, coords, grid_x, grid_y, u_1, v_1, u_2, v_2){

    # Replace each argument (a promise, which remembers the calling code) with
    # its plain value, so the frame holds only the grid
    type = type
    coords = coords
    grid_x = grid_x
    grid_y = grid_y
    u_1 = u_1
    v_1 = v_1
    u_2 = u_2
    v_2 = v_2

    # No source references in the returned function, so it stays small once
    # byte-compiled. It is not bound to a name here, so the call frame (its
    # environment) holds only the grid.
    utils::removeSource(switch(type,

      discontinuous = function(t, curPos){
        nearest_idx = which.min((curPos[1] - grid_x)^2 + (curPos[2] - grid_y)^2)
        matrix(c(u_1[nearest_idx], v_1[nearest_idx], u_2[nearest_idx], v_2[nearest_idx]), nrow = 2, byrow = FALSE)
      },

      continuous = function(t, curPos){
        interp = function(z){akima::bilinear(coords, coords, z, curPos[1], curPos[2])$z}
        matrix(c(interp(u_1), interp(v_1), interp(u_2), interp(v_2)), nrow = 2, byrow = FALSE)
      },

      differentiable = function(t, curPos){
        interp = function(z){akima::bicubic(coords, coords, z, curPos[1], curPos[2])$z}
        matrix(c(interp(u_1), interp(v_1), interp(u_2), interp(v_2)), nrow = 2, byrow = FALSE)
      }))

  }
  fieldFactory = utils::removeSource(fieldFactory)
  environment(fieldFactory) = baseenv()

  makeVectorFields = function(type, grid_1, grid_2){

    if(type == "discontinuous"){
      fieldFactory(type, coords = NULL, grid_x = grid_1$x, grid_y = grid_1$y,
                   u_1 = grid_1$u, v_1 = grid_1$v, u_2 = grid_2$u, v_2 = grid_2$v)
    } else{
      # Rows are x and columns are y, as expand.grid varies x fastest
      n = length(grid_coords)
      fieldFactory(type, coords = grid_coords, grid_x = NULL, grid_y = NULL,
                   u_1 = matrix(grid_1$u, nrow = n, ncol = n), v_1 = matrix(grid_1$v, nrow = n, ncol = n),
                   u_2 = matrix(grid_2$u, nrow = n, ncol = n), v_2 = matrix(grid_2$v, nrow = n, ncol = n))
    }

  }

  # ---- Settings --------------------------------------------------------------

  settings = expand.grid(error_level = names(error_sd), complexity = names(n_regions), continuity = continuity,
                         stringsAsFactors = FALSE)
  settings$error_sd = unname(error_sd[settings$error_level])
  settings$n_true_regions = unname(n_regions[settings$complexity])
  settings$name = stringr::str_c("error_", settings$error_level, ".complexity_", settings$complexity,
                                 ".continuity_", settings$continuity)

  # ---- Simulate one data set -------------------------------------------------

  simulateData = function(n_true_regions, sd, vectorFields){

    traj_splits = seq(0, t_end, length.out = n_true_regions + 1)

    # Redraw until the path stays on the grid, where the fields are defined
    for(attempt in seq_len(max_attempts)){

      true_traj = generateFullTrajectory(TimeSplits = traj_splits, numModels = n_models, sigma = sigma, k = k)
      start_pos = stats::runif(n_dim, min = start_range[1], max = start_range[2])

      true_data = RungeKutta(startTime = 0, startPos = start_pos, baseVectorFields = vectorFields,
                             TrajList = true_traj, TrajTimeSplits = traj_splits, endTime = t_end, t_step = t_step)

      if(all(abs(as.matrix(true_data[, -1])) <= 10)){
        break
      }

      if(attempt == max_attempts){
        stop("No path stayed inside the grid in ", max_attempts, " attempts.", call. = FALSE)
      }
    }

    sim_data = true_data
    for(d in seq_len(n_dim)){
      sim_data[, d + 1] = sim_data[, d + 1] + stats::rnorm(nrow(sim_data), mean = 0, sd = sd)
    }

    list(TrueTraj = list(Traj = true_traj, TimeSplits = traj_splits), TrueData = true_data, SimData = sim_data,
         Attempts = attempt)

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

  # Draw the replicate's grids, then simulate every data set, in a fixed order.
  set.seed(seed)
  grid_1 = createRandomGrid()
  grid_2 = createRandomGrid()

  vector_fields = lapply(stats::setNames(continuity, continuity), makeVectorFields, grid_1 = grid_1, grid_2 = grid_2)

  sim_data = lapply(seq_len(nrow(settings)), function(i){
    say("simulating data set ", i, "/", nrow(settings), ": ", settings$name[i])
    simulateData(n_true_regions = settings$n_true_regions[i], sd = settings$error_sd[i],
                 vectorFields = vector_fields[[settings$continuity[i]]])
  })
  names(sim_data) = settings$name

  results = list()
  fits = list()

  for(i in seq_len(nrow(settings))){

    fits[[settings$name[i]]] = list()

    for(variant in variants){

      fit = fitVariant(variant, sim_data[[i]], vector_fields[[settings$continuity[i]]], extra = list(path_attempts = sim_data[[i]]$Attempts))
      say("fit ", length(results) + 1, "/", nrow(settings)*length(variants), ": ", settings$name[i], ", ", variant,
          sprintf(" (%.1f s)", fit$metrics$time_total),
          if(!is.na(fit$metrics$error)) stringr::str_c(" FAILED: ", fit$metrics$error))

      results[[length(results) + 1]] = data.frame(seed = seed, settings[i, c("error_level", "error_sd", "complexity", "n_true_regions", "continuity")],
                                                  fit$metrics, row.names = NULL, check.names = FALSE)

      if(keep_fits){
        fits[[settings$name[i]]][[variant]] = list(EstPath = fit$EstPath, EstTraj = fit$EstTraj, Variance = fit$Variance)
      }
    }
  }

  # Attempts are reported in Results, not kept in Data
  sim_data = lapply(sim_data, function(s){s$Attempts = NULL; s})

  say("done in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), run_start, units = "mins"))), " min")

  out = list(Seed = seed,
             Settings = list(study = "Continuity2D", error_sd = error_sd, n_regions = n_regions, continuity = continuity,
                             variants = variants, B = B, patience = patience, t_end = t_end, t_step = t_step,
                             path_nsplits = path_nsplits, traj_nsplits = traj_nsplits, CL = CL,
                             sigma = sigma, k = k, start_range = start_range, max_attempts = max_attempts),
             Grids = list(grid_1, grid_2),
             VectorFields = vector_fields,
             Results = do.call(rbind, results),
             Data = sim_data)

  if(keep_fits){
    out$Fits = fits
  }

  out

}
