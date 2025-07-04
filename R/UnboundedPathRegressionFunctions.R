#' constructUnboundedCubicSplineModelMatrix_Path
#' 
#' @description
#' Constructing the model matrix for the path cubic spline regression. 
#' 
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions. 
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param split_labels a numeric or factor vector. Labels for each point in Data for which region of the spline the point falls into. 
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1. 
#'
#' @returns a numeric matrix. The regression model matrix X
#' @export
#'
constructUnboundedCubicSplineModelMatrix_Path = function(Data, time_splits, split_labels, V_Smooth = F){
  
  if(V_Smooth){
    
    a1 = ((Data$t - time_splits[1])^3)*(split_labels == 1) +
      ((time_splits[2] - time_splits[1])^3 + 3*(time_splits[2] - time_splits[1])*(Data$t - time_splits[1])*(Data$t - time_splits[2]))*(split_labels %in% c(2:(length(time_splits)-1)))
    b1 = (Data$t - time_splits[1])^2
    c1 = Data$t - time_splits[1]
    d1 = 1*(split_labels %in% c(1:(length(time_splits)-1)))
    
    cur_X = cbind(a1,b1,c1,d1)
    
    cur_X = cbind(cur_X, apply(matrix(2:(length(time_splits)-1)), MARGIN = 1, FUN = function(j){
      ((Data$t - time_splits[j])^3)*(split_labels == j) +
        ((time_splits[j+1] - time_splits[j])^3 + 3*(time_splits[j+1]-time_splits[j])*(Data$t - time_splits[j])*(Data$t - time_splits[j+1]))*(split_labels %in% c((j+1):(length(time_splits)-1)))*(j!=(length(time_splits)-1))
    }))
    
    cur_X
    
  } else{
    
    a1 = ((Data$t - time_splits[1])^3)*(split_labels == 1) +
      (3*(time_splits[2]-time_splits[1])^2*(Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^3)*(split_labels %in% c((2):(length(time_splits)-1)))
    b1 = ((Data$t - time_splits[1])^2)*(split_labels == 1) +
      (2*(time_splits[2]-time_splits[1])*(Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^2)*(split_labels %in% c((2):(length(time_splits)-1)))
    c1 = Data$t - time_splits[1]
    d1 = 1*(split_labels %in% c(1:(length(time_splits)-1)))
    
    cur_X = cbind(a1,b1,c1,d1)
    
    cur_X = cbind(cur_X, do.call("cbind", apply(matrix(2:(length(time_splits)-1)), MARGIN = 1, FUN = function(j){
      cur_a = ((Data$t - time_splits[j])^3)*(split_labels == j) +
        (3*(time_splits[j+1]-time_splits[j])^2*(Data$t - time_splits[j+1]) + (time_splits[j+1] - time_splits[j])^3)*(split_labels %in% c((j+1):(length(time_splits)-1)))*(j!=(length(time_splits)-1))
      cur_b = ((Data$t - time_splits[j])^2)*(split_labels == j) +
        (2*(time_splits[j+1]-time_splits[j])*(Data$t - time_splits[j+1]) + (time_splits[j+1] - time_splits[j])^2)*(split_labels %in% c((j+1):(length(time_splits)-1)))*(j!=(length(time_splits)-1))
      cbind(matrix(cur_a), matrix(cur_b))
    }, simplify = F)))
    
    cur_X
    
  }
  
}


#' Title
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions. 
#' @param dim a numeric scalar. The dimension of the data to run the regression. 
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param split_labels a numeric or factor vector. Labels for each point in Data for which region of the spline the point falls into. 
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#'
#' @returns either a list or NULL. NULL if the model matrix is singular. 
#' If the model matrix is non-singular the list contains the estimates of the coefficients, sigma2, and coefficient covariance matrix. 
#' @export
#'
runUnboundedCubicSplineRegression_Path = function(Data, dim, time_splits, split_labels, V_Smooth = F){

  cur_X = constructUnboundedCubicSplineModelMatrix_Path(Data, time_splits, split_labels, V_Smooth)
  cur_Y = Data[,dim+1]
  
  cur_reg = fastmatrix::ols.fit(x = cur_X, y = cur_Y)
  
  if(length(cur_reg$coefficients) != ncol(cur_X)){
    return(NULL)
  } else{
    
    cur_coef = unname(cur_reg$coefficients)
    cur_sigma2 = cur_reg$RSS/(nrow(cur_X) - ncol(cur_X))
    
    cur_cov = cur_reg$cov.unscaled*cur_sigma2
    
    list(cur_coef,cur_sigma2,cur_cov)
  }
  
  
}

