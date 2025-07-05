#' runEstTrajStep_Optimization
#'
#' @param EstPath an EtsimatedPath object.
#' @param baseVectorFields a baseVectorFields functions.
#' @param TimeSplits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param worst_split a numeric scalar. The region to be split.
#' @param cur_integral_list_traj a list of matrices. Each corresponds with a block of the current overall trajectory integral matrix.
#' @param cur_integral_list_pathtraj a list of matrices. Each corresponds with a block of the current overall PathTrajectory integral matrix.
#' @param n_dim a numeric scalar. The number of dimensions in the path.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#' @param Path_VSmooth a boolean. TRUE if the path was estimated using C2 smoothness. FALSE if C1.
#' @param Traj_VSmooth a boolean. TRUE if the trajectory should be estimated using C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedTrajectory object.
#' @export
#'
runEstTrajStep_Optimization = function(EstPath, baseVectorFields, TimeSplits, worst_split, cur_integral_list_traj, cur_integral_list_pathtraj, n_dim, n_models, Path_VSmooth, Traj_VSmooth){
  
  Path_TimeSplits = EstPath$TimeSplits
  
  free_path_coef = c(t(EstPath$Free_Parameters))
  
  free_path_cov = as.matrix(Matrix::bdiag(EstPath$Free_Parameter_Cov))
  
  Omega_p = constructUnboundedCubicSplineTransitionMatrix_Traj(Path_TimeSplits, n_dim, Path_VSmooth)
  Omega_t = constructUnboundedCubicSplineTransitionMatrix_Traj(TimeSplits, n_models, Traj_VSmooth)
  
  integral_list_traj = constructIntegralList_Traj(EstPath = EstPath, baseVectorFields = baseVectorFields, cur_integral_list = cur_integral_list_traj, time_splits = TimeSplits, worst_split = worst_split, n_models = n_models)
  Sigma_Traj = as.matrix(Matrix::bdiag(integral_list_traj))
  
  integral_list_pathtraj = constructIntegralList_PathTraj(EstPath = EstPath, baseVectorFields = baseVectorFields, cur_integral_list = cur_integral_list_pathtraj, worst_split = worst_split, Path_TimeSplits = Path_TimeSplits, Traj_TimeSplits = TimeSplits, n_dim = n_dim, n_models = n_models)
  Sigma_PathTraj = t(do.call(cbind, integral_list_pathtraj))
  
  Path_to_Traj_Matrix = Omega_t %*% pracma::pinv(t(Omega_t) %*% Sigma_Traj %*% Omega_t) %*% t(Omega_t) %*% Sigma_PathTraj %*% Omega_p
  
  traj_coef = t(Path_to_Traj_Matrix %*% free_path_coef)
  
  traj_coefs_list = lapply(unname(split(seq_len(ncol(traj_coef)), ceiling(seq_len(ncol(traj_coef)) / (4*n_models)))), function(cols) traj_coef[, cols, drop = FALSE])
  
  EstTraj = lapply(traj_coefs_list, FUN = function(coef){matrix(coef, ncol = 4, byrow = T)})
  
  full_traj_cov = Path_to_Traj_Matrix %*% free_path_cov %*% t(Path_to_Traj_Matrix)
  
  traj_cov_list = lapply(as.list(1:(length(TimeSplits) - 1)), FUN = function(x){full_traj_cov[1:(4*n_models) + (4*n_models)*(x-1), 1:(4*n_models) + (4*n_models)*(x-1)]})
  
  getSquaredDistancePerRegion = function(i){
    
    integrand = function(t){
      
      cur_pos = getEstPathPosition(t, EstPath)
      cur_vel = getEstPathVelocity(t, EstPath)
      
      est_vel = t(baseVectorFields(t, cur_pos) %*% t(getEstTrajValue(t = t, EstTraj = list(Traj = EstTraj, TimeSplits = TimeSplits))))
      
      ((cur_vel - est_vel) %*% t(cur_vel - est_vel))[1,1]
      
    }
    
    err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = TimeSplits[i], upper = TimeSplits[i+1], subdivisions = 10000)$value, 
                   error = function(e){NULL})
    
    if(is.null(err)){
      err = cubature::cubintegrate(integrand, lower = TimeSplits[i], upper = TimeSplits[i+1], method = 'pcubature')$integral
    }
    
    err
    
  }
  
  Error_v = apply(matrix(1:(length(TimeSplits)-1)), MARGIN = 1, FUN = getSquaredDistancePerRegion)
  Sigma2_list = replicate(length(TimeSplits)-1, list(sum(Error_v)/(TimeSplits[length(TimeSplits)]-TimeSplits[1])))
  
  path_split_labels_unique = cut(x = unique(EstPath$Error$t), breaks = TimeSplits, include.lowest = T, right = T, labels = 1:(length(TimeSplits)-1))
  
  pos_splits = unname(table(path_split_labels_unique) >= 8*n_models/n_dim)
  
  list(Traj = EstTraj, TimeSplits = TimeSplits, Error = Error_v, n_splits = length(TimeSplits)-1, Cov = traj_cov_list, Sigma2 = Sigma2_list, IntegralList_Traj = integral_list_traj, IntegralList_PathTraj = integral_list_pathtraj, PosSplits = pos_splits)
  
  
}



#' getEstimatedTrajectories_Optimization
#'
#' @param EstPath an EstimatedPath object.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#' @param baseVectorFields a baseVectorField object.
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#' @param Model_Selection_Type a string. The model criterion. 
#' Must pick from the following: BIC, BICOpt, BICPos, AIC, VarPos. 
#' @param Random_Path a boolean. TRUE if the path used to estimate the trajectories should be randomly generated. 
#' If TRUE, the path must've been estimated using the FullRegression or Projection method.
#'
#' @returns an EstimatedTrajectory object.
#' @export
#'
getEstimatedTrajectories_Optimization = function(EstPath, nsplits, baseVectorFields, V_Smooth = F, Model_Selection_Type = "BICOpt", Random_Path = F){
  
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  if(Random_Path){
    EstPath$Path = generateRandomPath(EstPath$Free_Parameters, free_param_cov = EstPath$Free_Parameter_Cov, TimeSplits = EstPath$TimeSplits, V_Smooth = V_Smooth)
  }
  
  n_dim = nrow(EstPath$Path[[1]])
  
  n_models = ncol(baseVectorFields(min(EstPath$TimeSplits),rep(0,n_dim)))
  
  time_splits = c(min(EstPath$TimeSplits), mean(c(min(EstPath$TimeSplits), max(EstPath$TimeSplits))), max(EstPath$TimeSplits))
  
  all_splits_est_traj = list()
  
  all_splits_est_traj[[1]] = runEstTrajStep_Optimization(EstPath = EstPath, baseVectorFields = baseVectorFields, TimeSplits = time_splits, worst_split = 1, cur_integral_list_traj = list(), cur_integral_list_pathtraj = list(), n_dim = n_dim, n_models = n_models, Path_VSmooth = V_Smooth, Traj_VSmooth = V_Smooth)
  
  svMisc::progress(1, max.value = nsplits, progress.bar = T, console = T)
  
  if(sum(all_splits_est_traj[[1]]$PosSplits) == 0){
    all_splits_est_traj
    
  } else{
    
    while(length(time_splits)-1 <= nsplits){
      
      split_labels_path = cut(x = EstPath$Error$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      median_path_error = unname(by(rowSums(EstPath$Error[,-1]), split_labels_path, stats::median))
      median_path_error[is.na(median_path_error)] = stats::median(median_path_error, na.rm = T)
      
      worst_split = which.max(((log(all_splits_est_traj[[length(all_splits_est_traj)]]$Error+1)) - median_path_error + max(median_path_error))*all_splits_est_traj[[length(all_splits_est_traj)]]$PosSplits)
      
      time_splits = c(time_splits[1:worst_split], mean(c(time_splits[worst_split], time_splits[worst_split+1])), time_splits[(worst_split+1):length(time_splits)])
      
      all_splits_est_traj[[length(all_splits_est_traj)+1]] = runEstTrajStep_Optimization(EstPath = EstPath, baseVectorFields = baseVectorFields, TimeSplits = time_splits, worst_split = worst_split, cur_integral_list_traj = all_splits_est_traj[[length(all_splits_est_traj)]]$IntegralList_Traj, cur_integral_list_pathtraj = all_splits_est_traj[[length(all_splits_est_traj)]]$IntegralList_PathTraj, n_dim = n_dim, n_models = n_models, Path_VSmooth = V_Smooth, Traj_VSmooth = V_Smooth)
      
      svMisc::progress(length(time_splits) - 2, max.value = nsplits, progress.bar = T, console = T)
      
      if(sum(all_splits_est_traj[[length(all_splits_est_traj)]]$PosSplits) == 0){
        break
      }
      
    }
  }
  
  getBestSplitsTraj(EstTrajList = all_splits_est_traj, EstPath = EstPath, show = T, type = Model_Selection_Type, V_Smooth = V_Smooth)
  
}


