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

  if(length(time_splits) == 2){
    # A single region: an ordinary cubic in (t - t_0).
    s = Data$t - time_splits[1]
    return(cbind(a1 = s^3, b1 = s^2, c1 = s, d1 = 1))
  }

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

#' constructUnboundedCubicSplineTransitionMatrix_Path
#' 
#' @description
#' Construct the matrix to transform the free parameters of the cubic spline into the full set of parameters.
#' 
#'
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#'
#' @returns a numeric matrix.
#' @export
#'
constructUnboundedCubicSplineTransitionMatrix_Path = function(time_splits, V_Smooth = T){
  
  # Each region's coefficients (a, b, c, d) follow from the previous region's
  # through the continuity conditions at the node between them, with
  # h = the previous region's length:
  #   d_j = a h^3 + b h^2 + c h + d    (position)
  #   c_j = 3 a h^2 + 2 b h + c        (velocity)
  #   b_j = 3 a h + b                  (acceleration, C2 only; free for C1)
  # The first region's four coefficients are free, then each later region adds
  # a_j (C2) or a_j and b_j (C1) as free parameters.
  
  n_regions = length(time_splits) - 1
  h = diff(time_splits)
  n_free = if(V_Smooth) 4 + (n_regions - 1) else 4 + 2*(n_regions - 1)
  
  trans_mat = matrix(0, 4*n_regions, n_free)
  trans_mat[1:4, 1:4] = diag(4)
  
  for(j in seq_len(n_regions)[-1]){
    
    prev = trans_mat[4*(j-2) + 1:4, , drop = FALSE]
    hj = h[j-1]
    rows = 4*(j-1) + 1:4
    
    if(V_Smooth){
      trans_mat[rows[1], 4 + (j-1)] = 1
      trans_mat[rows[2], ] = 3*hj*prev[1, ] + prev[2, ]
    } else{
      trans_mat[rows[1], 4 + 2*(j-2) + 1] = 1
      trans_mat[rows[2], 4 + 2*(j-2) + 2] = 1
    }
    
    trans_mat[rows[3], ] = 3*hj^2*prev[1, ] + 2*hj*prev[2, ] + prev[3, ]
    trans_mat[rows[4], ] = hj^3*prev[1, ] + hj^2*prev[2, ] + hj*prev[3, ] + prev[4, ]
    
  }
  
  trans_mat
}


#' runUnboundedCubicSplineRegression_Path
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

