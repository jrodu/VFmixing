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
  
  if(V_Smooth){
    
    trans_mat = matrix(c(1,rep(0,length(time_splits)+1),
                         0,1,rep(0,length(time_splits)),
                         0,0,1,rep(0,length(time_splits)-1),
                         0,0,0,1,rep(0,length(time_splits)-2)), byrow = T, nrow = 4)
    
    creation_Func = function(i){
      if(i%%4 == 1){
        
        cur_row = c(rep(0,4 + (i-5)/4),1,rep(0,length(time_splits) - 3 - (i-5)/4))
        
      } else if(i%%4 == 2){
        
        sub_time_splits = time_splits[2:((i+2)/4)]
        
        v_init = c(3*diff(time_splits)[1],1,0,0)
        
        v1 = 3*diff(sub_time_splits)
        
        cur_row = c(v_init, v1, rep(0,length(time_splits) - (i+2)/4))
        
      } else if(i%%4 == 3){
        
        sub_time_splits = time_splits[2:((i+1)/4)]
        
        v_init = c(3*diff(time_splits)[1]^2 + 6*diff(time_splits)[1]*(time_splits[(i+1)/4] - time_splits[2]),
                   2*(time_splits[(i+1)/4] - time_splits[1]),
                   1,
                   0)
        
        v1 = 3*diff(sub_time_splits)^2 + 6*diff(sub_time_splits)*(time_splits[(i+1)/4] - sub_time_splits[-1])
        
        cur_row = c(v_init, v1, rep(0,length(time_splits) - (i+1)/4))
        
      }else if(i%%4 == 0){
        
        sub_time_splits = time_splits[2:(i/4)]
        
        v_init = c(diff(time_splits)[1]^3 + 3*diff(time_splits)[1]*(time_splits[i/4] - time_splits[2])*(time_splits[i/4] - time_splits[1]), 
                   (time_splits[i/4]-time_splits[1])^2,
                   time_splits[i/4]-time_splits[1],
                   1)
        
        v1 = diff(sub_time_splits)^3 + 3*diff(sub_time_splits)^2*(time_splits[i/4] - sub_time_splits[-1]) + 3*diff(sub_time_splits)*(time_splits[i/4] - sub_time_splits[-1])^2
        
        cur_row = c(v_init, v1,rep(0,length(time_splits) - i/4))
      }
      
      cur_row
    }
    
    trans_mat = rbind(trans_mat, t(apply(matrix(5:(4*(length(time_splits)-1))), MARGIN = 1, FUN = creation_Func, simplify = T)))
    
  } else{
    
    trans_mat = matrix(c(1,rep(0,2*length(time_splits)-1),
                         0,1,rep(0,2*length(time_splits)-2),
                         0,0,1,rep(0,2*length(time_splits)-3),
                         0,0,0,1,rep(0,2*length(time_splits)-4)), byrow = T, nrow = 4)
    
    creation_Func = function(i){
      
      if(i%%4 == 1){
        cur_row = c(rep(0,(i+3)/2),1,rep(0,2*length(time_splits) - ((i+3)/2) - 1))
      } else if(i%%4 == 2){
        cur_row = c(rep(0, (i+4)/2),1,rep(0,2*length(time_splits) - ((i+4)/2) - 1))
      } else if(i%%4 == 3){
        
        sub_time_splits = time_splits[2:((i+1)/4)]
        
        v_init = c(3*diff(time_splits)[1]^2,2*diff(time_splits)[1],1,0)
        
        v1 = 3*diff(sub_time_splits)^2
        v2 = 2*diff(sub_time_splits)
        
        cur_row = c(v_init, c(rbind(v1,v2)),rep(0,2*length(time_splits) - 4 - 2*(length(sub_time_splits)-1)))
        
      }else if(i%%4 == 0){
        
        sub_time_splits = time_splits[2:(i/4)]
        
        v_init = c(diff(time_splits)[1]^3 + 3*diff(time_splits)[1]^2*(time_splits[i/4] - time_splits[2]), 
                   diff(time_splits)[1]^2 + 2*diff(time_splits)[1]*(time_splits[i/4] - time_splits[2]),
                   time_splits[i/4]-time_splits[1],
                   1)
        
        v1 = diff(sub_time_splits)^3 + 3*diff(sub_time_splits)^2*(time_splits[i/4] - sub_time_splits[-1])
        v2 = diff(sub_time_splits)^2 + 2*diff(sub_time_splits)*(time_splits[i/4] - sub_time_splits[-1])
        
        cur_row = c(v_init, c(rbind(v1,v2)),rep(0,2*length(time_splits) - 4 - 2*(length(sub_time_splits)-1)))
      }
      
      cur_row
      
    }
    
    trans_mat = rbind(trans_mat, t(apply(matrix(5:(4*(length(time_splits)-1))), MARGIN = 1, FUN = creation_Func, simplify = T)))
    
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

