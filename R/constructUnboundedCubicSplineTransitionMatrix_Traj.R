#' constructUnboundedCubicSplineTransitionMatrix_Traj
#'
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#' @param V_Smooth a boolean. TRUE if the trajectory should be estimated using C2 smoothness. FALSE if C1.
#'
#' @returns a numeric matrix. Transforms the free trajectory coefficients into the full coefficients.
#' @export
#'
constructUnboundedCubicSplineTransitionMatrix_Traj = function(time_splits, n_models, V_Smooth = T){
  
  Path_TransMat = constructUnboundedCubicSplineTransitionMatrix_Path(time_splits, V_Smooth)
  
  Traj_TransMat = c()
  
  for(i in 1:(length(time_splits)-1)){
    
    Traj_TransMat = rbind(Traj_TransMat, as.matrix(Matrix::bdiag(replicate(n_models, Path_TransMat[1:4 + (i-1)*4,], simplify = F))))
    
  }
  
  Traj_TransMat
  
}
