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
#' cubicPercentError(EstCubic = ExEstTraj_Optimization$Traj, 
#'                   TrueCubic = ExEstTraj_ByFours$Traj, 
#'                   EstTimeSplits = ExEstTraj_Optimization$TimeSplits, 
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
#' cubicAltPercentError(EstCubic = ExEstTraj_Optimization$Traj, 
#'                      TrueCubic = ExTraj$Traj, 
#'                      EstTimeSplits = ExEstTraj_Optimization$TimeSplits, 
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
#' cubicAvgSquareDist(EstCubic = ExEstTraj_Optimization$Traj, 
#'                    TrueCubic = ExTraj$Traj, 
#'                    EstTimeSplits = ExEstTraj_Optimization$TimeSplits, 
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
#' @param EstTraj an EstimatedTrajectory object. Must have been estimated using the Optimzation method.
#' @param TrueTraj a list of numeric matrices. The list of coefficient matrices of a true trajectory.
#' @param TrueTrajTimeSplits a numeric vector. The nodes of the true trajectory in addition with the starting and ending time.
#' @param CL a real number between 0 and 1. The confidence level of the intervals to check the coverage rate.
#'
#' @returns a numeric vector. Each entry corresponds to the coverage rate for a given model
#' @export
#' 
#' @examples
#' coverageRateTraj(EstTraj = ExEstTraj_Optimization, TrueTraj = ExTraj$Traj, 
#'                  TrueTrajTimeSplits = ExTraj$TimeSplits, CL = 0.95)
#'   
#'
coverageRateTraj = function(EstTraj, TrueTraj, TrueTrajTimeSplits, CL = 0.95){
  
  TimeSplits = EstTraj$TimeSplits
  n_models = nrow(EstTraj$Traj[[1]])
  Cov = EstTraj$Cov
  
  t_seq = seq(min(TimeSplits), max(TimeSplits), by = 0.01)
  n_points = length(t_seq)
  
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
      
      cur_y_var = (t(trans_mat) %*% Cov[[cur_est_traj_index]][1:4 + (j-1)*4, 1:4 + (j-1)*4] %*% trans_mat)[1,1] + (EstTraj$Sigma2[[cur_est_traj_index]])
      upper = est_traj[j] + stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      lower = est_traj[j] - stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      
      inInterval = c(inInterval, (true_traj[j] > lower) & (true_traj[j] < upper))
      
    }
    
    
    inInterval
    
  }
  
  colMeans(matrix((apply(matrix(1:length(t_seq)), MARGIN = 1, FUN = isInInterval))))
  
  
  
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
#' avgMagError(baseVectorFields = baseVectorFields, Traj = ExEstTraj_Optimization, 
#'             Path = ExEstPath_Projection, upper = 5, lower = 0)
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