#' generateRandomPath
#' 
#'
#' @param free_param_mean a numeric matrix. 
#' Each column corresponds to the means of the free parameters for a dimension's estimated path cubic spline.
#' @param free_param_cov a list of numeric matrices. 
#' Each element corresponds to the covariance matrix of the free parameters for a dimension's estimated path cubic spline
#' @param TimeSplits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#'
#' @returns a list of matrices. Each element corresponds to the cubic spline coefficients for each region.
#' @export
#' 
#' @examples
#' generateRandomPath(free_param_mean = ExEstPath_FullRegression$Free_Parameters, 
#'                    free_param_cov = ExEstPath_FullRegression$Free_Parameter_Cov,
#'                    TimeSplits = ExEstPath_FullRegression$TimeSplits, V_Smooth = FALSE)
#' 
generateRandomPath = function(free_param_mean, free_param_cov, TimeSplits, V_Smooth){
  
  n_dim = nrow(free_param_mean)
  
  V = constructUnboundedCubicSplineTransitionMatrix_Path(TimeSplits, V_Smooth)
  
  random_FP = c()
  
  for(i in 1:n_dim){
    
    random_FP = rbind(random_FP, MASS::mvrnorm(n = 1, mu = free_param_mean[i,], Sigma = free_param_cov[[i]]))
    
  }
  
  
  random_FullCoef = t(V %*% t(random_FP))
  
  lapply(unname(split(seq_len(ncol(random_FullCoef)), ceiling(seq_len(ncol(random_FullCoef)) / 4))), function(cols) random_FullCoef[, cols, drop = FALSE])
  
}