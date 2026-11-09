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
#' **Refining the trajectory (`traj_refine`).** With low positional error and
#' complex trajectories, most of the remaining bias comes from trajectory nodes
#' that sit away from where the true trajectories change. The path can only
#' support so many regions, but the trajectory can be split further: each
#' selected trajectory region is split into `traj_refine` equal regions and refit
#' (see the `refine` argument of [getEstimatedTrajectories_OneShot()]). This is a
#' heuristic that trades wider intervals for less bias. In the bootstrap, every
#' replicate's trajectory is refined the same way.
#'
#' **Bootstrap (`bootstrap`).** OneShot's covariance treats the selected nodes as
#' fixed. The nodes are chosen from the data, so they vary from data set to data
#' set, and that variability is missing from it. By the law of total variance,
#' \deqn{Var(\hat\beta) = E[Var(\hat\beta \mid nodes)] + Var(E[\hat\beta \mid nodes]).}
#' The first term is OneShot's covariance (the within-node part). The second
#' (the selection part) is estimated with a residual bootstrap:
#' 1. The positions of the fitted path at the data times are taken as the truth,
#'    with the centred residuals as the noise. The residuals are scaled by
#'    \eqn{\sqrt{n/(n - p)}}, with \eqn{p} the path's free parameters per
#'    dimension, so the resampled errors have the same variance as the path's
#'    own estimate, \eqn{RSS/(n - p)}: least-squares residuals are smaller than
#'    the errors.
#' 2. For each of `B` replicates, the residuals are resampled with replacement
#'    and added to the fitted positions, and the whole method (both node
#'    searches) is rerun, giving that replicate's path and trajectory nodes.
#' 3. The replicate's nodes are then fit, without any noise, to the fitted
#'    positions. This is the conditional mean \eqn{m_b = E[\hat\beta \mid nodes]}
#'    for those nodes, a C1 cubic spline on the replicate's own trajectory nodes.
#' 4. The selection part at each time is the sample variance (divisor
#'    \eqn{B - 1}) of the conditional means there,
#'    \eqn{V_{sel}(t) = \frac{1}{B - 1}\sum_b (m_b(t) - \bar m(t))^2}.
#'    The conditional means are kept as they are, with no projection onto the
#'    original nodes.
#'
#' The pointwise intervals are then
#' \eqn{\hat\beta(t) \pm z \sqrt{s_w^2(t) + V_{sel}(t)}}, with \eqn{s_w^2(t)} the
#' within-node variance from OneShot's covariance. The returned trajectory keeps
#' OneShot's covariance as `Cov` and carries the conditional means as
#' `SelectionMeans`; [visualizeEstTraj()], [coverageRateTraj()] and
#' [avgCIWidthTraj()] add the selection part whenever `SelectionMeans` is
#' present. Set it to `NULL` for intervals from OneShot's covariance alone.
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
#' @param traj_refine a positive whole number. Each selected trajectory region is split into this many equal
#' regions and refit. The default, 1, keeps the selected regions. See Details.
#' @param bootstrap a boolean. If TRUE, add the node-selection variance to the trajectory covariance. See Details. Default FALSE.
#' @param B a positive integer. The number of bootstrap replicates.
#' @param quad_nodes a positive integer. Passed to [getEstimatedTrajectories_OneShot()].
#' @param t_grid_size a positive number. The time grid for the pointwise variances in `Variance$Pointwise`.
#' @param quiet a boolean. If TRUE, the progress output of the estimation functions is suppressed.
#'
#' @returns a list with:
#' * `EstPath`: the EstimatedPath object.
#' * `EstTraj`: the EstimatedTrajectory object, with OneShot's covariance as `Cov`. If `bootstrap`, it also has
#'   `SelectionMeans`: the conditional means, one per successful replicate (each a list with `Traj` and `TimeSplits`).
#' * `Variance`: if `bootstrap`, a list with `Within` (OneShot's covariance per region) and `Pointwise`, a data
#'   frame with the time, trajectory, and the within, selection and total variances of the trajectory value with
#'   the selection share of the total. NULL otherwise.
#' * `Bootstrap`: if `bootstrap`, a list with `B`, `Successful` (the number of replicates that fit), `Nodes`
#'   (each replicate's path and trajectory nodes) and `ResidualScale` (the residuals' scale factor). NULL otherwise.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fit = getEstimatedTrajectories_Full(ExData, baseVectorFields, traj_nsplits = 50,
#'                                     bootstrap = TRUE, B = 10, quiet = TRUE)
#' head(fit$Variance$Pointwise)
#' # Coverage of the true trajectories with the total variance
#' coverageRateTraj(fit$EstTraj, ExTraj$Traj, ExTraj$TimeSplits)
#' }
#'
getEstimatedTrajectories_Full = function(Data, baseVectorFields, path_nsplits = NULL, traj_nsplits = 100,
                                         path_method = "FullFours", V_Smooth = F, Model_Selection_Type = "BIC",
                                         patience = 10, undersmooth = T, traj_refine = 1, bootstrap = F, B = 20,
                                         quad_nodes = 20, t_grid_size = 0.01, quiet = F){

  if(!(path_method %in% c("FullFours", "FullRegression", "Projection"))){
    stop("path_method must be \"FullFours\", \"FullRegression\" or \"Projection\": OneShot needs a path with a covariance structure.",
         call. = FALSE)
  }

  if(!(Model_Selection_Type %in% c("BIC", "AIC"))){
    stop("Model_Selection_Type must be \"BIC\" or \"AIC\".", call. = FALSE)
  }

  checkPatience(patience)
  checkRefine(traj_refine)

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
                                     patience = patience, undersmooth = undersmooth, refine = traj_refine)
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

  # Least-squares residuals are smaller than the errors: rescale them so the
  # resampled errors have variance RSS/(n - p), with p the path's free
  # parameters per dimension, matching the path's own variance estimate.
  n_obs = nrow(resid)
  p_path = splineFreeParams(length(EstPath$Path), inferPathVSmooth(EstPath, default = V_Smooth))
  residual_scale = if(n_obs > p_path) sqrt(n_obs/(n_obs - p_path)) else 1
  resid = resid*residual_scale

  means = vector("list", B)
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

    if(!is.null(cond_traj)){
      # Kept as it is, on the replicate's own nodes: no projection
      means[[b]] = list(Traj = cond_traj$Traj, TimeSplits = cond_traj$TimeSplits)
    }
  }

  ok = !vapply(means, is.null, logical(1))
  n_ok = sum(ok)

  if(n_ok < 2){
    stop("Fewer than two bootstrap replicates could be fit.", call. = FALSE)
  }

  # OneShot's covariance stays the within part; the conditional means give the
  # selection part pointwise, and the interval functions add the two.
  EstTraj$SelectionMeans = means[ok]

  list(EstPath = EstPath,
       EstTraj = EstTraj,
       Variance = list(Within = EstTraj$Cov, Pointwise = pointwiseTrajVariance(EstTraj, t_grid_size)),
       Bootstrap = list(B = B, Successful = n_ok, Nodes = nodes[ok], ResidualScale = residual_scale))

}

#' Refine the Trajectories of a Full-Method Fit
#'
#' @description
#' Splits every trajectory region of a [getEstimatedTrajectories_Full()] fit
#' into `refine` equal regions and refits, reusing the fit's node selections
#' instead of searching again. With a bootstrap, every replicate's conditional
#' mean is refit on its own refined trajectory nodes too, so the result keeps
#' the node-selection variance.
#'
#' @details
#' Refinement happens after node selection, so every node search (on the data
#' and in each bootstrap replicate) would select the same nodes with or without
#' it. This gives the same result as rerunning [getEstimatedTrajectories_Full()]
#' with `traj_refine` multiplied by `refine` on the same bootstrap replicates,
#' at the cost of fixed-node refits only. Comparing a fit with its refined
#' version is also a paired comparison: both use the same replicates.
#'
#' @param fit a list returned by [getEstimatedTrajectories_Full()].
#' @param Data the data the fit was made from.
#' @param baseVectorFields the baseVectorField function the fit was made with.
#' @param refine a positive whole number. Each trajectory region is split into this many equal regions.
#' @param V_Smooth,quad_nodes,t_grid_size the values the fit was made with (the defaults match those of
#' [getEstimatedTrajectories_Full()]).
#'
#' @returns a list like that of [getEstimatedTrajectories_Full()], with the refined trajectories.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fit = getEstimatedTrajectories_Full(ExData, baseVectorFields, traj_nsplits = 50,
#'                                     bootstrap = TRUE, B = 10, quiet = TRUE)
#' refined = refineEstimatedTrajectories_Full(fit, ExData, baseVectorFields, refine = 2)
#' length(fit$EstTraj$Traj)
#' length(refined$EstTraj$Traj)
#' }
#'
refineEstimatedTrajectories_Full = function(fit, Data, baseVectorFields, refine = 2, V_Smooth = F, quad_nodes = 20,
                                            t_grid_size = 0.01){

  checkRefine(refine)

  if(refine == 1){
    return(fit)
  }

  Data = data.frame(Data)
  names(Data) = c("t", stringr::str_c("X", seq_len(ncol(Data) - 1)))

  EstTraj = oneShotFixedNodes(fit$EstPath, refineTimeSplits(fit$EstTraj$TimeSplits, refine), baseVectorFields, V_Smooth, quad_nodes)

  if(is.null(fit$Bootstrap)){
    return(list(EstPath = fit$EstPath, EstTraj = EstTraj, Variance = NULL, Bootstrap = NULL))
  }

  # Each replicate's conditional mean on its own refined trajectory nodes, with
  # its path nodes fit to the noise-free fitted positions as before
  Data_fit = Data
  Data_fit[, -1] = getEstPathPosition(Data$t, fit$EstPath)

  nodes = lapply(fit$Bootstrap$Nodes, function(n){list(Path = n$Path, Traj = refineTimeSplits(n$Traj, refine))})

  EstTraj$SelectionMeans = lapply(nodes, function(n){
    EstPath_m = runEstPathStep_FullRegression(Data_fit, n$Path, V_Smooth = V_Smooth)
    cond_traj = oneShotFixedNodes(EstPath_m, n$Traj, baseVectorFields, V_Smooth, quad_nodes)
    list(Traj = cond_traj$Traj, TimeSplits = cond_traj$TimeSplits)
  })

  Bootstrap = fit$Bootstrap
  Bootstrap$Nodes = nodes

  list(EstPath = fit$EstPath,
       EstTraj = EstTraj,
       Variance = list(Within = EstTraj$Cov, Pointwise = pointwiseTrajVariance(EstTraj, t_grid_size)),
       Bootstrap = Bootstrap)

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

# The within-node variance of each trajectory's value at times t, from
# OneShot's covariance (one row per time, one column per trajectory).
trajWithinVariance = function(EstTraj, t){

  TimeSplits = EstTraj$TimeSplits
  n_models = nrow(EstTraj$Traj[[1]])
  idx = pmin(pmax(findInterval(t, TimeSplits, left.open = TRUE), 1), length(TimeSplits) - 1)
  u = t - TimeSplits[idx]
  B = cbind(u^3, u^2, u, 1)

  variance = matrix(0, length(t), n_models)
  for(r in unique(idx)){
    sel = idx == r
    Bs = B[sel, , drop = FALSE]
    for(m in seq_len(n_models)){
      block = EstTraj$Cov[[r]][1:4 + (m - 1)*4, 1:4 + (m - 1)*4]
      variance[sel, m] = rowSums((Bs %*% block) * Bs)
    }
  }
  pmax(variance, 0)

}

# The node-selection variance of each trajectory's value at times t: the sample
# variance (divisor B - 1) of the bootstrap conditional means at each time.
# Zero if the trajectory has no bootstrap (one row per time, one column per trajectory).
trajSelectionVariance = function(EstTraj, t){

  n_models = nrow(EstTraj$Traj[[1]])

  if(is.null(EstTraj$SelectionMeans)){
    return(matrix(0, length(t), n_models))
  }

  values = vapply(EstTraj$SelectionMeans, function(m){evaluateTrajAt(t, m)}, matrix(0, length(t), n_models))
  values = array(values, c(length(t), n_models, length(EstTraj$SelectionMeans)))
  apply(values, c(1, 2), stats::var)

}

# Within, selection and total variances of each trajectory's value on a time grid.
pointwiseTrajVariance = function(EstTraj, t_grid_size){

  TimeSplits = EstTraj$TimeSplits
  t_grid = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)
  within = trajWithinVariance(EstTraj, t_grid)
  selection = trajSelectionVariance(EstTraj, t_grid)

  rows = lapply(seq_len(ncol(within)), function(m){
    total = within[, m] + selection[, m]
    data.frame(t = t_grid, trajectory = m, var_within = within[, m], var_selection = selection[, m],
               var_total = total, selection_share = ifelse(total > 0, selection[, m]/total, 0))
  })

  do.call(rbind, rows)

}
