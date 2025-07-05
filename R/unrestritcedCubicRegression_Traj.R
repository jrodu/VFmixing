#' constructUnrestrictedCubicRegressionModelMatrix_Traj
#'
#' @param cur_PathData a numeric data frame. The path data to use to estimate the cubic trajectories. 
#' @param start_time a numeric scalar. The centering time of the cubic function.
#' @param n_dim a numeric scalar. The number of dimensions in the path.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#'
#' @returns a numeric matrix. The linear regression model matrix.
#' @export
#'
constructUnrestrictedCubicRegressionModelMatrix_Traj = function(cur_PathData, start_time, n_dim, n_models){
  
  cur_X = c()
  
  for(i in 1:n_dim){
    
    tmp_X = c()
    
    for(j in 1:n_models){
      
      a1 = (cur_PathData$t - start_time)^3
      b1 = (cur_PathData$t - start_time)^2
      c1 = (cur_PathData$t - start_time)
      d1 = 1
      
      tmp_X = cbind(tmp_X, cur_PathData[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(a1,b1,c1,d1))
      
    }
    
    cur_X = rbind(cur_X, tmp_X)
    
  }
  
  cur_X
  
}


#' runUnRestrictedCubicRegression_Traj
#'
#' @param cur_PathData a numeric data frame. The path data to use to estimate the cubic trajectories. 
#' @param start_time a numeric scalar. The centering time of the cubic function.
#' @param n_dim a numeric scalar. The number of dimensions in the path.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#'
#' @returns a list. The elements are as follows: the estimated coefficients, the sigma2 estimate, the model matrix.
#' @export
#'
runUnRestrictedCubicRegression_Traj = function(cur_PathData, start_time, n_dim, n_models){
  
  cur_X = constructUnrestrictedCubicRegressionModelMatrix_Traj(cur_PathData = cur_PathData, start_time = start_time, n_dim = n_dim, n_models = n_models)  
  
  cur_Y = c()
  
  for(i in 1:n_dim){
    
    cur_Y = c(cur_Y, cur_PathData[,1+n_dim+i])
    
  }
  
  cur_coef = t(backsolve(qr.R(qr(cur_X)), t(qr.Q(qr(cur_X))) %*% cur_Y))
  
  cur_sigma2 = mean(((cur_X %*% t(cur_coef)) - cur_Y)^2)
  
  list(cur_coef,cur_sigma2, cur_X)
  
  
}

