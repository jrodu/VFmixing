#' Full Method: Estimated Path and Trajectories with Optional Bootstrap
#'
#' @description
#' Runs the whole method on the data: a path is estimated with a covariance
#' structure (FullFours by default), then the trajectories are estimated from it
#' with [getEstimatedTrajectories_OneShot()]. Both node searches use the "2x BIC"
#' undersmoothing rule by default. Optionally, a bootstrap adds the uncertainty
#' from the node selection to OneShot's covariance through the law of total
#' variance.
#'
#' @details
#' **Undersmoothing (`undersmooth`).** Each node search returns the candidate
#' with twice as many regions as the one its criterion selects. This reduces the
#' bias of the path's velocity and of the trajectory spline, so the trajectory
#' intervals are centred closer to the truth. See [getEstimatedPaths_FullFours()].
#'
#' **Bootstrap (`bootstrap`).** OneShot's covariance treats the selected nodes as
#' fixed. The nodes are chosen from the data, so they vary from data set to data
#' set, and that variability is missing from it. By the law of total variance,
#' \deqn{Var(\hat\beta) = E[Var(\hat\beta \mid nodes)] + Var(E[\hat\beta \mid nodes]).}
#' The first term is OneShot's covariance (the within-node part). The second
#' (the selection part) is estimated with a residual bootstrap:
#' 1. The positions of the fitted path at the data times are taken as the truth,
#'    with the centred residuals as the noise.
#' 2. For each of `B` replicates, the residuals are resampled with replacement
#'    and added to the fitted positions, and the whole method (both node
#'    searches) is rerun, giving that replicate's path and trajectory nodes.
#' 3. The replicate's nodes are then fit, without any noise, to the fitted
#'    positions. This is the conditional mean \eqn{E[\hat\beta \mid nodes]} for
#'    those nodes.
#' 4. The selection part is the variance of those conditional means across
#'    replicates. They are projected onto the original trajectory's regions
#'    (a cubic per region, by least squares), so the selection covariance has the
#'    same form as OneShot's and the two can be added.
#'
#' The returned trajectory's `Cov` is the total (within + selection), so
#' [visualizeEstTraj()], [coverageRateTraj()] and [avgCIWidthTraj()] use the
#' total uncertainty. The separate pieces are returned in `Variance`.
#'
#' Because the conditional means are fit without noise, the selection part
#' settles with few replicates: 20 is usually plenty. This part still conditions
#' on the base vector fields evaluated along the estimated path, so the error
#' from the path position in those evaluations is not included.
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param baseVectorFields a baseVectorField function.
#' @param path_nsplits a positive integer, or NULL for a quarter of the number of data points.
#' The maximum number of splits for the path search.
#' @param traj_nsplits a positive integer. The maximum number of splits for the trajectory search.
#' @param path_method a string. The path method: "FullFours", "FullRegression" or "Projection".
#' @param V_Smooth a boolean. TRUE for C2 smoothness of the path and trajectory splines, FALSE for C1.
#' @param Model_Selection_Type a string. The criterion for both node searches, "BIC" or "AIC".
#' @param patience NULL or a positive integer. Passed to both node searches.
#' @param undersmooth a boolean. If TRUE (the default), both node searches use twice the selected number of regions. See Details.
#' @param bootstrap a boolean. If TRUE, add the node-selection variance to the trajectory covariance. See Details. Default FALSE.
#' @param B a positive integer. The number of bootstrap replicates.
#' @param quad_nodes a positive integer. Passed to [getEstimatedTrajectories_OneShot()].
#' @param t_grid_size a positive number. The time grid for the pointwise variances in `Variance$Pointwise`.
#' @param quiet a boolean. If TRUE, the progress output of the estimation functions is suppressed.
#'
#' @returns a list with:
#' * `EstPath`: the EstimatedPath object.
#' * `EstTraj`: the EstimatedTrajectory object. Its `Cov` is the total covariance (OneShot's covariance plus the
#'   selection part if `bootstrap`, otherwise OneShot's covariance alone).
#' * `Variance`: if `bootstrap`, a list with `Within` (OneShot's covariance per region), `Selection` (the
#'   selection covariance per region), `Total` (their sum) and `Pointwise`, a data frame with the time, trajectory,
#'   and the within, selection and total variances of the trajectory value with the selection share of the total.
#'   NULL otherwise.
#' * `Bootstrap`: if `bootstrap`, a list with `B`, `Successful` (the number of replicates that fit) and `Nodes`
#'   (each replicate's path and trajectory nodes). NULL otherwise.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fit = getEstimatedTrajectories_Full(ExData, baseVectorFields, traj_nsplits = 50,
#'                                     bootstrap = TRUE, B = 10, quiet = TRUE)
#' head(fit$Variance$Pointwise)
#' # Coverage of the true trajectories with the total covariance
#' coverageRateTraj(fit$EstTraj, ExTraj$Traj, ExTraj$TimeSplits)
#' }
#'
getEstimatedTrajectories_Full = function(Data, baseVectorFields, path_nsplits = NULL, traj_nsplits = 100,
                                         path_method = "FullFours", V_Smooth = F, Model_Selection_Type = "BIC",
                                         patience = 10, undersmooth = T, bootstrap = F, B = 20,
                                         quad_nodes = 20, t_grid_size = 0.01, quiet = F){

  if(!(path_method %in% c("FullFours", "FullRegression", "Projection"))){
    stop("path_method must be \"FullFours\", \"FullRegression\" or \"Projection\": OneShot needs a path with a covariance structure.",
         call. = FALSE)
  }

  if(!(Model_Selection_Type %in% c("BIC", "AIC"))){
    stop("Model_Selection_Type must be \"BIC\" or \"AIC\".", call. = FALSE)
  }

  checkPatience(patience)

  if(bootstrap && (length(B) != 1 || !is.numeric(B) || B < 2 || B != round(B))){
    stop("B must be a whole number of at least 2.", call. = FALSE)
  }

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c("t", stringr::str_c("X", 1:n_dim))

  if(is.null(path_nsplits)){
    path_nsplits = floor(nrow(Data)/4)
  }

  runQuietly = function(expr, silent){
    if(silent){
      res = NULL
      utils::capture.output(res <- expr)
      res
    } else{
      expr
    }
  }

  fitPath = function(d){
    switch(path_method,
           FullFours = getEstimatedPaths_FullFours(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                   patience = patience, undersmooth = undersmooth),
           FullRegression = getEstimatedPaths_FullRegression(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                             patience = patience, undersmooth = undersmooth),
           Projection = getEstimatedPaths_Projection(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                     patience = patience, undersmooth = undersmooth))
  }

  fitTraj = function(EstPath){
    getEstimatedTrajectories_OneShot(EstPath, nsplits = traj_nsplits, baseVectorFields = baseVectorFields, V_Smooth = V_Smooth,
                                     Model_Selection_Type = Model_Selection_Type, quad_nodes = quad_nodes,
                                     patience = patience, undersmooth = undersmooth)
  }

  # ---- The fit on the data ---------------------------------------------------

  EstPath = runQuietly(fitPath(Data), quiet)
  EstTraj = runQuietly(fitTraj(EstPath), quiet)

  if(!bootstrap){
    return(list(EstPath = EstPath, EstTraj = EstTraj, Variance = NULL, Bootstrap = NULL))
  }

  # ---- The selection part of the variance ------------------------------------

  fitted = getEstPathPosition(Data$t, EstPath)
  Data_fit = Data
  Data_fit[, -1] = fitted
  resid = as.matrix(Data[, -1]) - fitted
  resid = sweep(resid, 2, colMeans(resid))

  TimeSplits = EstTraj$TimeSplits
  n_regions = length(TimeSplits) - 1
  n_models = nrow(EstTraj$Traj[[1]])

  # Times on each original trajectory region, to project the conditional means onto its cubics
  region_t = lapply(seq_len(n_regions), function(r){seq(TimeSplits[r], TimeSplits[r + 1], length.out = 25)})

  coef_draws = lapply(seq_len(n_regions), function(r){matrix(NA_real_, B, 4*n_models)})
  nodes = vector("list", B)

  for(b in seq_len(B)){

    if(!quiet){
      message("Bootstrap replicate ", b, "/", B)
    }

    Data_b = Data_fit
    Data_b[, -1] = fitted + resid[sample(nrow(resid), replace = TRUE), , drop = FALSE]

    cond_traj = tryCatch(runQuietly({
      EstPath_b = fitPath(Data_b)
      EstTraj_b = fitTraj(EstPath_b)
      nodes[[b]] = list(Path = EstPath_b$TimeSplits, Traj = EstTraj_b$TimeSplits)

      # The conditional mean for this replicate's nodes: the same nodes fit to the noise-free fitted positions
      EstPath_m = runEstPathStep_FullRegression(Data_fit, EstPath_b$TimeSplits, V_Smooth = V_Smooth)
      oneShotFixedNodes(EstPath_m, EstTraj_b$TimeSplits, baseVectorFields, V_Smooth, quad_nodes)
    }, TRUE), error = function(e){NULL})

    if(is.null(cond_traj)){
      next
    }

    for(r in seq_len(n_regions)){
      u = region_t[[r]] - TimeSplits[r]
      X = cbind(u^3, u^2, u, 1)
      y = evaluateTrajAt(region_t[[r]], cond_traj)
      coef_draws[[r]][b, ] = c(qr.solve(X, y))
    }
  }

  ok = !is.na(coef_draws[[1]][, 1])
  n_ok = sum(ok)

  if(n_ok < 2){
    stop("Fewer than two bootstrap replicates could be fit.", call. = FALSE)
  }

  Cov_Within = EstTraj$Cov
  Cov_Selection = lapply(coef_draws, function(draws){stats::cov(draws[ok, , drop = FALSE])})
  Cov_Total = mapply(`+`, Cov_Within, Cov_Selection, SIMPLIFY = FALSE)

  EstTraj$Cov = Cov_Total

  list(EstPath = EstPath,
       EstTraj = EstTraj,
       Variance = list(Within = Cov_Within, Selection = Cov_Selection, Total = Cov_Total,
                       Pointwise = pointwiseTrajVariance(TimeSplits, Cov_Within, Cov_Selection, n_models, t_grid_size)),
       Bootstrap = list(B = B, Successful = n_ok, Nodes = nodes[ok]))

}

# OneShot trajectory with fixed nodes (no node search) on a given path.
oneShotFixedNodes = function(EstPath, TimeSplits, baseVectorFields, V_Smooth = F, quad_nodes = 20){

  n_dim = nrow(EstPath$Path[[1]])
  n_models = ncol(baseVectorFields(min(EstPath$TimeSplits), rep(0, n_dim)))
  gl = gaussLegendre(quad_nodes)

  Omega_p = constructUnboundedCubicSplineTransitionMatrix_Traj(EstPath$TimeSplits, n_dim, inferPathVSmooth(EstPath, default = V_Smooth))

  caches = lapply(seq_len(length(TimeSplits) - 1), function(r){
    oneShotRegionCache(TimeSplits[r], TimeSplits[r + 1], EstPath, baseVectorFields, n_dim, n_models, gl)
  })

  runEstTrajStep_OneShot(EstPath, TimeSplits, caches, Omega_p, n_dim, n_models, V_Smooth)

}

# The trajectory values (one column per model) at times t.
evaluateTrajAt = function(t, EstTraj){

  ts = EstTraj$TimeSplits
  idx = pmin(pmax(findInterval(t, ts, left.open = TRUE), 1), length(ts) - 1)
  u = t - ts[idx]
  n_models = nrow(EstTraj$Traj[[1]])

  out = matrix(NA_real_, length(t), n_models)
  for(i in seq_along(t)){
    out[i, ] = EstTraj$Traj[[idx[i]]] %*% c(u[i]^3, u[i]^2, u[i], 1)
  }
  out

}

# Within, selection and total variances of each trajectory's value on a time grid.
pointwiseTrajVariance = function(TimeSplits, Cov_Within, Cov_Selection, n_models, t_grid_size){

  t_grid = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)
  idx = pmin(pmax(findInterval(t_grid, TimeSplits, left.open = TRUE), 1), length(TimeSplits) - 1)

  rows = lapply(seq_len(n_models), function(m){
    k = 1:4 + (m - 1)*4
    v = t(sapply(seq_along(t_grid), function(i){
      u = t_grid[i] - TimeSplits[idx[i]]
      b = c(u^3, u^2, u, 1)
      c(drop(t(b) %*% Cov_Within[[idx[i]]][k, k] %*% b), drop(t(b) %*% Cov_Selection[[idx[i]]][k, k] %*% b))
    }))
    data.frame(t = t_grid, trajectory = m, var_within = v[, 1], var_selection = v[, 2],
               var_total = v[, 1] + v[, 2], selection_share = v[, 2]/(v[, 1] + v[, 2]))
  })

  do.call(rbind, rows)

}
