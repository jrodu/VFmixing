#' cubicPercentError
#' 
#' @description
#' Get the percent error between 2 cubic spline functions.
#'  
#'
#' @param EstCubic a list of numeric matrices. The list of coefficient matrices of the estimated cubic spline.
#' @param TrueCubic a list of numeric matrices. The list of coefficient matrices of a true cubic spline.
#' @param EstTimeSplits a numeric vector. The nodes of the estimated cubic spline in addition with the starting and ending time.
#' @param TrueTimeSplits a numeric vector. The nodes of the true cubic spline in addition with the starting and ending time.
#'
#' @returns a numeric vector. The first is the percent error and the second is the error associated with the percent error.
#' @export
#' 
#' @examples
#' cubicPercentError(EstCubic = ExEstTraj_OneShot$Traj, 
#'                   TrueCubic = ExEstTraj_ByFours$Traj, 
#'                   EstTimeSplits = ExEstTraj_OneShot$TimeSplits, 
#'                   TrueTimeSplits = ExEstTraj_ByFours$TimeSplits)
#'   
cubicPercentError = function(EstCubic, TrueCubic, EstTimeSplits, TrueTimeSplits){
  
  integrand = function(t){
    
    cur_est_traj_index = min(c(max(c(sum((EstTimeSplits - t)<0),1)), length(EstTimeSplits)-1))
    cur_true_traj_index = min(c(max(c(sum((TrueTimeSplits - t)<0),1)), length(TrueTimeSplits)-1))
    curEstTraj = EstCubic[[cur_est_traj_index]]
    curTrueTraj = TrueCubic[[cur_true_traj_index]]
    
    curEst = evaluateCubic(t = t, Cubic = curEstTraj, startTime = EstTimeSplits[cur_est_traj_index])
    curTrue = evaluateCubic(t = t, Cubic = curTrueTraj, startTime = TrueTimeSplits[cur_true_traj_index])
    
    mean(abs((curEst - curTrue) / (curTrue)))
  }
  
  lower_bound = min(c(EstTimeSplits, TrueTimeSplits))
  upper_bound = max(c(EstTimeSplits, TrueTimeSplits))
  
  err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = lower_bound, upper = upper_bound, subdivisions = 10000, rel.tol = 0.00001), 
                 error = function(e){NULL})
  
  if(is.null(err)){
    err = cubature::cubintegrate(integrand, lower = lower_bound, upper = upper_bound, method = 'pcubature')
    
    return_v = c(Value = 100*err$integral / (upper_bound - lower_bound), Error = 100*err$error / (upper_bound - lower_bound))
    
  } else{
    
    return_v = c(Value = 100*err$value / (upper_bound - lower_bound), Error = 100*err$abs.error / (upper_bound - lower_bound))
  }
  
  return_v
  
}

#' cubicAltPercentError
#' 
#' @description
#' Another cubic percent error-like metric.
#' 
#'
#' @param EstCubic a list of numeric matrices. The list of coefficient matrices of the estimated cubic spline.
#' @param TrueCubic a list of numeric matrices. The list of coefficient matrices of a true cubic spline.
#' @param EstTimeSplits a numeric vector. The nodes of the estimated cubic spline in addition with the starting and ending time.
#' @param TrueTimeSplits a numeric vector. The nodes of the true cubic spline in addition with the starting and ending time.
#'
#' @returns a numeric scalar. 
#' @export
#' 
#' @examples
#' cubicAltPercentError(EstCubic = ExEstTraj_OneShot$Traj, 
#'                      TrueCubic = ExTraj$Traj, 
#'                      EstTimeSplits = ExEstTraj_OneShot$TimeSplits, 
#'                      TrueTimeSplits = ExTraj$TimeSplits)
#' 
cubicAltPercentError = function(EstCubic, TrueCubic, EstTimeSplits, TrueTimeSplits){
  
  altPE = c()
  
  for(i in 1:nrow(EstCubic[[1]])){
    
    integrand = function(t){
      cur_est_traj_index = min(c(max(c(sum((EstTimeSplits - t)<0),1)),length(EstTimeSplits)-1))
      cur_true_traj_index = min(c(max(c(sum((TrueTimeSplits - t)<0),1)),length(TrueTimeSplits)-1))
      curEstTraj = EstCubic[[cur_est_traj_index]][i,]
      curTrueTraj = TrueCubic[[cur_true_traj_index]][i,]
      
      curEst = evaluateCubic(t = t, Cubic = curEstTraj, startTime = EstTimeSplits[cur_est_traj_index])
      curTrue = evaluateCubic(t = t, Cubic = curTrueTraj, startTime = TrueTimeSplits[cur_true_traj_index])
      
      abs(curEst - curTrue)
    }
    
    lower_bound = min(c(EstTimeSplits, TrueTimeSplits))
    upper_bound = max(c(EstTimeSplits, TrueTimeSplits))
    
    absdiff = stats::integrate(f = Vectorize(integrand), lower = lower_bound, upper = upper_bound, subdivisions = 10000, rel.tol = 0.001)$value
    
    integrand = function(t){
      cur_true_traj_index = min(c(max(c(sum((TrueTimeSplits - t)<0),1)),length(TrueTimeSplits)-1))
      curTrueTraj = TrueCubic[[cur_true_traj_index]][i,]
      
      curTrue = evaluateCubic(t = t, Cubic = curTrueTraj, startTime = TrueTimeSplits[cur_true_traj_index])
      
      abs(curTrue)
    }
    
    truth = stats::integrate(f = Vectorize(integrand), lower = lower_bound, upper = upper_bound, subdivisions = 10000, rel.tol = 0.001)$value
    
    altPE = c(altPE, absdiff/truth)
    
  }
  
  
  mean(altPE)
  
}

#' cubicAvgSquareDist
#' 
#' @description
#' Finds the average squared distance between 2 cubic spline functions.
#' 
#'
#' @param EstCubic a list of numeric matrices. The list of coefficient matrices of the estimated cubic spline.
#' @param TrueCubic a list of numeric matrices. The list of coefficient matrices of a true cubic spline.
#' @param EstTimeSplits a numeric vector. The nodes of the estimated cubic spline in addition with the starting and ending time.
#' @param TrueTimeSplits a numeric vector. The nodes of the true cubic spline in addition with the starting and ending time.
#' @param lower_bound a number scalar. The lower bound of the intergal.
#' @param upper_bound a number scalar. The upper bound of the intergal.
#'
#' @returns a numeric scalar.
#' @export
#' 
#' @examples
#' cubicAvgSquareDist(EstCubic = ExEstTraj_OneShot$Traj, 
#'                    TrueCubic = ExTraj$Traj, 
#'                    EstTimeSplits = ExEstTraj_OneShot$TimeSplits, 
#'                    TrueTimeSplits = ExTraj$TimeSplits)
#'   
cubicAvgSquareDist = function(EstCubic, TrueCubic, EstTimeSplits, TrueTimeSplits, lower_bound, upper_bound){
  
  integrand = function(t){
    cur_est_traj_index = min(c(max(c(sum((EstTimeSplits - t)<0),1)),length(EstTimeSplits)-1))
    cur_true_traj_index = min(c(max(c(sum((TrueTimeSplits - t)<0),1)),length(TrueTimeSplits)-1))
    curEstTraj = EstCubic[[cur_est_traj_index]]
    curTrueTraj = TrueCubic[[cur_true_traj_index]]
    
    curEst = evaluateCubic(t = t, Cubic = curEstTraj, startTime = EstTimeSplits[cur_est_traj_index])
    curTrue = evaluateCubic(t = t, Cubic = curTrueTraj, startTime = TrueTimeSplits[cur_true_traj_index])
    
    mean((curEst - curTrue)^2)
  }
  
  if(missing(lower_bound) | missing(upper_bound)){
    
    lower_bound = min(c(EstTimeSplits, TrueTimeSplits))
    upper_bound = max(c(EstTimeSplits, TrueTimeSplits))
    
  }
  
  
  
  err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = lower_bound, upper = upper_bound, subdivisions = 10000)$value / 
                   (upper_bound - lower_bound), 
                 error = function(e){NULL})
  
  if(is.null(err)){
    err = cubature::cubintegrate(integrand, lower = lower_bound, upper = upper_bound, method = 'pcubature')$integral / (upper_bound - lower_bound) 
    
  }
  
  err
}

#' coverageRatePath
#'
#' @param EstPath an EstimatedPath object. Must have been estimated using either the FullRegression or Projection methods.
#' @param Data a numeric data frame.
#' @param CL a real number between 0 and 1. The confidence level of the intervals to check the coverage rate.
#'
#' @returns a numeric vector. Each entry corresponds to the coverage rate for a given dimension.
#' @export
#' 
#' @examples
#' coverageRatePath(ExEstPath_Projection, ExData, 0.95)
#' 
coverageRatePath = function(EstPath, Data, CL = 0.95){
  
  TimeSplits = EstPath$TimeSplits
  n_dim = nrow(EstPath$Path[[1]])
  Cov = EstPath$Cov
  
  isInInterval = function(i){
    
    t = Data$t[i]
    
    truePos = Data[i,-1]
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    est_pos = evaluateCubic(t, EstPath$Path[[cur_path_index]], startTime = TimeSplits[cur_path_index])
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
    
    inInterval = c()
    
    for(j in 1:n_dim){
      
      cur_y_var = (t(trans_mat) %*% Cov[[cur_path_index]][[j]] %*% trans_mat)[1,1] + (EstPath$Sigma2[[cur_path_index]][j])
      upper = est_pos[j] + stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      lower = est_pos[j] - stats::qnorm((1+CL)/2)*sqrt(cur_y_var) 
      
      inInterval = c(inInterval, (truePos[j] > lower) & (truePos[j] < upper))
      
    }
    
    inInterval
    
  }
  
  
  colMeans(t(apply(matrix(1:length(Data$t)), MARGIN = 1, FUN = isInInterval)))
  
  
  
}

#' coverageRateTraj
#'
#' @param EstTraj an EstimatedTrajectory object. Must have been estimated using the OneShot method.
#' @param TrueTraj a list of numeric matrices. The list of coefficient matrices of a true trajectory.
#' @param TrueTrajTimeSplits a numeric vector. The nodes of the true trajectory in addition with the starting and ending time.
#' @param CL a real number between 0 and 1. The confidence level of the intervals to check the coverage rate.
#' The intervals are confidence intervals from the trajectory coefficient covariance.
#'
#' @returns a numeric vector. Each entry corresponds to the coverage rate for a given model
#' @export
#' 
#' @examples
#' coverageRateTraj(EstTraj = ExEstTraj_OneShot, TrueTraj = ExTraj$Traj, 
#'                  TrueTrajTimeSplits = ExTraj$TimeSplits, CL = 0.95)
#'   
#'
coverageRateTraj = function(EstTraj, TrueTraj, TrueTrajTimeSplits, CL = 0.95){
  
  TimeSplits = EstTraj$TimeSplits
  n_models = nrow(EstTraj$Traj[[1]])
  Cov = EstTraj$Cov
  
  t_seq = seq(min(TimeSplits), max(TimeSplits), by = 0.01)
  n_points = length(t_seq)

  # The node-selection variance, if the trajectory has bootstrap conditional
  # means (getEstimatedTrajectories_Full()); zero otherwise
  selection_var = trajSelectionVariance(EstTraj, t_seq)

  isInInterval = function(i){
    
    t = t_seq[i]
    
    cur_est_traj_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_true_traj_index = min(c(max(c(sum((TrueTrajTimeSplits - t)<0),1)), length(TrueTrajTimeSplits)-1))
    est_traj = evaluateCubic(t, EstTraj$Traj[[cur_est_traj_index]], startTime = TimeSplits[cur_est_traj_index])
    true_traj = evaluateCubic(t, TrueTraj[[cur_true_traj_index]], startTime = TrueTrajTimeSplits[cur_true_traj_index])
    
    cur_difft = t - TimeSplits[cur_est_traj_index]
    trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
    
    inInterval = c()
    
    for(j in 1:n_models){
      
      # Confidence interval: the trajectories are not observed, so no noise term is added
      cur_y_var = (t(trans_mat) %*% Cov[[cur_est_traj_index]][1:4 + (j-1)*4, 1:4 + (j-1)*4] %*% trans_mat)[1,1] + selection_var[i, j]
      upper = est_traj[j] + stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      lower = est_traj[j] - stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      
      inInterval = c(inInterval, (true_traj[j] > lower) & (true_traj[j] < upper))
      
    }
    
    
    inInterval
    
  }
  
  # apply() returns one column per time point and one row per model.
  inInterval_mat = matrix(apply(matrix(1:length(t_seq)), MARGIN = 1, FUN = isInInterval), nrow = n_models)
  
  rowMeans(inInterval_mat)
  
  
  
}


#' avgMagError
#' 
#' @description
#' Finds the average magnitude of the error velocity vectors. 
#' 
#'
#' @param baseVectorFields a baseVectorFields function.
#' @param Traj an EstimatedTrajectory object.
#' @param Path an EstimatedPath object.
#' @param upper a numeric scalar. The upper bound for the integral.
#' @param lower a numeric scalar. The lower bound for the integral.
#'
#' @returns a numeric scalar.
#' @export
#' 
#' @examples
#' avgMagError(baseVectorFields = baseVectorFields, Traj = ExEstTraj_OneShot, 
#'             Path = ExEstPath_FullFours, upper = 5, lower = 0)
#'   
avgMagError = function(baseVectorFields, Traj, Path, upper, lower){
  
  integrand = function(t){
    
    curPos = getEstPathPosition(t, Path)
    curVel = getEstPathVelocity(t, Path)
    WeightedVF = TrajWeightedBaseVectorFields(t, curPos, baseVectorFields, Traj$Traj, Traj$TimeSplits)
    
    sum((WeightedVF - curVel)^2)/sum(curVel^2)
    
  }
  
  err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = lower, upper = upper, subdivisions = 10000, rel.tol = 0.00001)$value, 
                 error = function(e){NULL})
  
  if(is.null(err)){
    
    err = cubature::cubintegrate(integrand, lower = lower, upper = upper, method = 'pcubature')$integral
    
  }
  
  err/(upper-lower)
  
  
}

#' avgCIWidthPath
#'
#' @description
#' The average width of an estimated path's confidence (or prediction)
#' intervals over time, for each dimension.
#'
#' @details
#' At time \eqn{t} the interval for dimension \eqn{d} has width
#' \eqn{2 z \sqrt{b(t)^T \Sigma_d b(t)}}, where \eqn{b(t)} is the cubic basis,
#' \eqn{\Sigma_d} the coefficient covariance and \eqn{z} the normal quantile for
#' `CL`. With `P_Int = TRUE` the observation noise \eqn{\sigma^2_d} is added
#' under the square root, which gives the width of the prediction interval used
#' by [coverageRatePath()] and [visualizeEstPath()]. The widths are averaged
#' over an evenly spaced time grid, so this is a time-average across the whole
#' path.
#'
#' @param EstPath an EstimatedPath object with a covariance structure (FullFours, FullRegression or Projection).
#' @param CL a real number between 0 and 1. The confidence level of the intervals.
#' @param P_Int a boolean. FALSE (the default) for confidence intervals on the path,
#' TRUE for prediction intervals (path uncertainty plus observation noise).
#' @param t_grid_size a positive numeric scalar. The spacing of the time grid the widths are averaged over.
#'
#' @returns a named numeric vector. The average interval width for each dimension.
#' @export
#'
#' @examples
#' avgCIWidthPath(ExEstPath_FullFours, CL = 0.95)
#' avgCIWidthPath(ExEstPath_FullFours, CL = 0.95, P_Int = TRUE)
#'
avgCIWidthPath = function(EstPath, CL = 0.95, P_Int = F, t_grid_size = 0.01){

  checkIntervalArgs(CL, t_grid_size)

  if(is.null(EstPath$Cov) || is.null(EstPath$Sigma2)){
    stop("EstPath has no covariance structure. Use a FullFours, FullRegression or Projection path.", call. = FALSE)
  }

  TimeSplits = EstPath$TimeSplits
  n_dim = nrow(EstPath$Path[[1]])
  t_seq = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)

  path = evaluatePathAt(t_seq, EstPath)
  B = cbind(path$s^3, path$s^2, path$s, 1)
  variance = matrix(0, length(t_seq), n_dim)

  for(r in unique(path$idx)){
    sel = path$idx == r
    Bs = B[sel, , drop = FALSE]
    for(d in seq_len(n_dim)){
      variance[sel, d] = rowSums((Bs %*% EstPath$Cov[[r]][[d]]) * Bs) + P_Int*EstPath$Sigma2[[r]][d]
    }
  }

  z = stats::qnorm((1 + CL)/2)
  width = colMeans(2*z*sqrt(pmax(variance, 0)))
  names(width) = stringr::str_c("X", seq_len(n_dim))
  width

}


#' avgCIWidthTraj
#'
#' @description
#' The average width of an estimated trajectory's confidence intervals over
#' time, for each trajectory.
#'
#' @details
#' At time \eqn{t} the interval for trajectory \eqn{m} has width
#' \eqn{2 z \sqrt{b(t)^T \Sigma_m b(t)}}, where \eqn{b(t)} is the cubic basis,
#' \eqn{\Sigma_m} that trajectory's block of the coefficient covariance and
#' \eqn{z} the normal quantile for `CL`. These are the intervals drawn by
#' [visualizeEstTraj()] and checked by [coverageRateTraj()]. The widths are
#' averaged over an evenly spaced time grid, so this is a time-average across
#' the whole trajectory.
#'
#' There is no prediction-interval option: the trajectories are never observed
#' directly, and the OneShot \eqn{\sigma^2} is the variance of the velocity
#' residuals, not of a trajectory value.
#'
#' @param EstTraj an EstimatedTrajectory object with a covariance structure (OneShot).
#' @param CL a real number between 0 and 1. The confidence level of the intervals.
#' @param t_grid_size a positive numeric scalar. The spacing of the time grid the widths are averaged over.
#'
#' @returns a named numeric vector. The average interval width for each trajectory.
#' @export
#'
#' @examples
#' avgCIWidthTraj(ExEstTraj_OneShot, CL = 0.95)
#'
avgCIWidthTraj = function(EstTraj, CL = 0.95, t_grid_size = 0.01){

  checkIntervalArgs(CL, t_grid_size)

  if(is.null(EstTraj$Cov)){
    stop("EstTraj has no covariance structure. Use a OneShot trajectory.", call. = FALSE)
  }

  TimeSplits = EstTraj$TimeSplits
  n_models = nrow(EstTraj$Traj[[1]])
  t_seq = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)

  # OneShot's variance, plus the node-selection variance if the trajectory has
  # bootstrap conditional means (getEstimatedTrajectories_Full())
  variance = trajWithinVariance(EstTraj, t_seq) + trajSelectionVariance(EstTraj, t_seq)

  z = stats::qnorm((1 + CL)/2)
  width = colMeans(2*z*sqrt(pmax(variance, 0)))
  names(width) = stringr::str_c("Trajectory", seq_len(n_models))
  width

}


# Shared argument checks for the interval-width functions.
checkIntervalArgs = function(CL, t_grid_size){

  if(length(CL) != 1 || !is.numeric(CL) || CL <= 0 || CL >= 1){
    stop("CL must be a number between 0 and 1.", call. = FALSE)
  }

  if(length(t_grid_size) != 1 || !is.numeric(t_grid_size) || t_grid_size <= 0){
    stop("t_grid_size must be a positive number.", call. = FALSE)
  }

  invisible(NULL)
}
