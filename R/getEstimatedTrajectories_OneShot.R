#' runEstTrajStep_OneShot
#'
#' @description
#' One candidate in the OneShot trajectory method: the trajectory spline for a
#' given set of nodes. Called by [getEstimatedTrajectories_OneShot()].
#'
#' @param EstPath an EstimatedPath object.
#' @param TimeSplits a numeric vector. The nodes of the trajectory spline with the starting and ending times.
#' @param RegionCaches a list with one element per trajectory region, as built internally by
#' [getEstimatedTrajectories_OneShot()]. Each holds the quadrature nodes, the base vector fields and
#' path velocity at those nodes, and the region's blocks of the integral matrices.
#' @param Omega_p a numeric matrix. The transition matrix from the path free parameters to the full path coefficients.
#' @param n_dim a numeric scalar. The number of dimensions in the path.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#' @param Traj_VSmooth a boolean. TRUE if the trajectory should be estimated using C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedTrajectory object.
#' @export
#'
runEstTrajStep_OneShot = function(EstPath, TimeSplits, RegionCaches, Omega_p, n_dim, n_models, Traj_VSmooth){

  n_regions = length(TimeSplits) - 1

  free_path_coef = c(t(EstPath$Free_Parameters))
  free_path_cov = as.matrix(Matrix::bdiag(EstPath$Free_Parameter_Cov))

  Omega_t = trajTransitionMatrix(TimeSplits, n_models, Traj_VSmooth)

  integral_list_traj = lapply(RegionCaches, `[[`, "G")
  Sigma_Traj = blockDiagDense(integral_list_traj)
  Sigma_PathTraj = do.call(rbind, lapply(RegionCaches, `[[`, "S"))

  # Least-squares map from the path coefficients to the trajectory coefficients:
  # Omega_t (Omega_t' Sigma_Traj Omega_t)^+ Omega_t' Sigma_PathTraj Omega_p
  Path_to_Traj_Matrix = Omega_t %*% psdSolve(crossprod(Omega_t, Sigma_Traj) %*% Omega_t,
                                             crossprod(Omega_t, Sigma_PathTraj %*% Omega_p))

  traj_coef = c(Path_to_Traj_Matrix %*% free_path_coef)

  EstTraj = lapply(seq_len(n_regions), function(r){
    matrix(traj_coef[(r - 1)*4*n_models + 1:(4*n_models)], ncol = 4, byrow = T)
  })

  full_traj_cov = Path_to_Traj_Matrix %*% tcrossprod(free_path_cov, Path_to_Traj_Matrix)

  traj_cov_list = lapply(seq_len(n_regions), function(x){
    full_traj_cov[1:(4*n_models) + (4*n_models)*(x-1), 1:(4*n_models) + (4*n_models)*(x-1)]
  })

  # Integrated squared velocity error in each region
  Error_v = vapply(seq_len(n_regions), function(r){
    oneShotRegionError(RegionCaches[[r]], EstTraj[[r]], n_dim, n_models)
  }, numeric(1))

  Sigma2_list = replicate(n_regions, list(sum(Error_v)/(TimeSplits[length(TimeSplits)]-TimeSplits[1])))

  path_split_labels_unique = cut(x = unique(EstPath$Error$t), breaks = TimeSplits, include.lowest = T, right = T, labels = 1:n_regions)

  pos_splits = unname(table(path_split_labels_unique) >= 8*n_models/n_dim)

  list(Traj = EstTraj, TimeSplits = TimeSplits, Error = Error_v, n_splits = n_regions, Cov = traj_cov_list, Sigma2 = Sigma2_list,
       IntegralList_Traj = integral_list_traj,
       IntegralList_PathTraj = lapply(RegionCaches, function(cache){t(cache$S)}),
       PosSplits = pos_splits)

}



#' getEstimatedTrajectories_OneShot
#'
#' @description
#' Estimates the trajectory (weight) functions by transforming the path
#' coefficients into trajectory coefficients, which gives a covariance
#' structure for the trajectories.
#'
#' @details
#' The integrals this method needs are computed with Gauss-Legendre quadrature:
#' each trajectory region is split at the path nodes, and the path and
#' `baseVectorFields` are evaluated at `quad_nodes` points on each piece. Every
#' integral is then formed from those values, and a region's values are reused
#' for every later candidate. With smooth vector fields the default of 20 nodes
#' is accurate to around machine precision. Increase `quad_nodes` if
#' `baseVectorFields` changes abruptly along the path.
#'
#' Candidates are produced by repeatedly splitting a region, up to `nsplits`
#' splits. By default every candidate is fit and the best is chosen with
#' `Model_Selection_Type`. Setting `patience` stops the search once that many
#' candidates in a row have failed to improve on the best criterion so far.
#' This is much faster when the best model has few regions, but it can miss a
#' better model that would only appear after a long run without improvement.
#'
#' @param EstPath an EstimatedPath object.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#' @param baseVectorFields a baseVectorField object.
#' @param V_Smooth a boolean. TRUE if the trajectory spline should have C2 smoothness. FALSE if C1.
#' The path's smoothness is detected from EstPath, so it does not need to match.
#' @param Model_Selection_Type a string. The model criterion.
#' Must pick from the following: BIC, BICOS, AIC, VarPos. See [getBestSplitsTraj()].
#' @param Random_Path a boolean. TRUE if the path used to estimate the trajectories should be randomly generated.
#' If TRUE, the path must've been estimated using the FullRegression, Projection or FullFours method.
#' @param quad_nodes a positive integer. The number of Gauss-Legendre nodes per piece of each region. See Details.
#' @param patience NULL or a positive integer. If NULL (the default), every candidate up to `nsplits` is fit.
#' Otherwise, stop after this many consecutive candidates without improving the criterion. See Details.
#'
#' @returns an EstimatedTrajectory object.
#' @export
#'
#' @examples
#' \donttest{
#' ExEstTraj_OneShot = getEstimatedTrajectories_OneShot(
#'      EstPath = ExEstPath_FullFours,
#'      nsplits = 100, baseVectorFields = baseVectorFields,
#'      V_Smooth = FALSE, Model_Selection_Type = "BIC", Random_Path = FALSE)
#' }
#'
#' # Stop once 10 candidates in a row fail to improve the BIC
#' EstTraj = getEstimatedTrajectories_OneShot(
#'      EstPath = ExEstPath_FullFours,
#'      nsplits = 100, baseVectorFields = baseVectorFields, patience = 10)
#'
getEstimatedTrajectories_OneShot = function(EstPath, nsplits, baseVectorFields, V_Smooth = F, Model_Selection_Type = "BIC", Random_Path = F, quad_nodes = 20, patience = NULL){

  Model_Selection_Type = checkTrajSelectionType(Model_Selection_Type)

  if(length(quad_nodes) != 1 || !is.numeric(quad_nodes) || quad_nodes < 1 || quad_nodes != round(quad_nodes)){
    stop("quad_nodes must be a positive whole number.", call. = FALSE)
  }

  checkPatience(patience)

  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)

  # The path may have been fit with different smoothness than the trajectory.
  Path_V_Smooth = inferPathVSmooth(EstPath, default = V_Smooth)

  if(Random_Path){
    EstPath$Path = generateRandomPath(EstPath$Free_Parameters, free_param_cov = EstPath$Free_Parameter_Cov, TimeSplits = EstPath$TimeSplits, V_Smooth = Path_V_Smooth)
  }

  n_dim = nrow(EstPath$Path[[1]])

  n_models = ncol(baseVectorFields(min(EstPath$TimeSplits),rep(0,n_dim)))

  gl = gaussLegendre(quad_nodes)

  # The path side never changes between candidates.
  Omega_p = constructUnboundedCubicSplineTransitionMatrix_Traj(EstPath$TimeSplits, n_dim, Path_V_Smooth)

  newRegion = function(lower, upper){
    oneShotRegionCache(lower, upper, EstPath, baseVectorFields, n_dim, n_models, gl)
  }

  fitCandidate = function(time_splits, caches){
    runEstTrajStep_OneShot(EstPath = EstPath, TimeSplits = time_splits, RegionCaches = caches, Omega_p = Omega_p, n_dim = n_dim, n_models = n_models, Traj_VSmooth = V_Smooth)
  }

  candidateCriterion = function(cand){
    trajSelectionCriterion(list(cand), EstPath, Model_Selection_Type, V_Smooth)
  }

  time_splits = c(min(EstPath$TimeSplits), mean(c(min(EstPath$TimeSplits), max(EstPath$TimeSplits))), max(EstPath$TimeSplits))

  caches = list(newRegion(time_splits[1], time_splits[2]), newRegion(time_splits[2], time_splits[3]))

  all_splits_est_traj = list()

  all_splits_est_traj[[1]] = fitCandidate(time_splits, caches)

  if(!is.null(patience)){
    best_crit = candidateCriterion(all_splits_est_traj[[1]])
    since_best = 0
  }

  svMisc::progress(1, max.value = nsplits, progress.bar = T, console = T)

  if(sum(all_splits_est_traj[[1]]$PosSplits) > 0){

    while(length(time_splits)-1 <= nsplits){

      last = all_splits_est_traj[[length(all_splits_est_traj)]]

      split_labels_path = cut(x = EstPath$Error$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      median_path_error = unname(by(rowSums(EstPath$Error[,-1]), split_labels_path, stats::median))
      median_path_error[is.na(median_path_error)] = stats::median(median_path_error, na.rm = T)

      worst_split = which.max(((log(last$Error+1)) - median_path_error + max(median_path_error))*last$PosSplits)

      mid_time = mean(c(time_splits[worst_split], time_splits[worst_split+1]))

      # Only the two new regions need new quadrature values.
      caches = append(caches[-worst_split],
                      list(newRegion(time_splits[worst_split], mid_time), newRegion(mid_time, time_splits[worst_split+1])),
                      after = worst_split - 1)

      time_splits = c(time_splits[1:worst_split], mid_time, time_splits[(worst_split+1):length(time_splits)])

      all_splits_est_traj[[length(all_splits_est_traj)+1]] = fitCandidate(time_splits, caches)

      svMisc::progress(length(time_splits) - 2, max.value = nsplits, progress.bar = T, console = T)

      if(sum(all_splits_est_traj[[length(all_splits_est_traj)]]$PosSplits) == 0){
        break
      }

      if(!is.null(patience)){
        cur_crit = candidateCriterion(all_splits_est_traj[[length(all_splits_est_traj)]])
        if(cur_crit < best_crit){
          best_crit = cur_crit
          since_best = 0
        } else{
          since_best = since_best + 1
          if(since_best >= patience){
            break
          }
        }
      }

    }
  }

  getBestSplitsTraj(EstTrajList = all_splits_est_traj, EstPath = EstPath, show = F, type = Model_Selection_Type, V_Smooth = V_Smooth)

}
