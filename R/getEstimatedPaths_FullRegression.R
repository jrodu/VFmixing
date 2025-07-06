#' runEstPathStep_FullRegression
#' 
#' @description
#' One iteration in the FullRegression method for estimated paths.
#' 
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions. 
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedPath object.
#' @export
#'
runEstPathStep_FullRegression = function(Data, time_splits, V_Smooth = T){
  
  est_path = vector(mode = "list", length = length(time_splits)-1)
  
  n_dim = ncol(Data)-1
  
  sigma2_list = list()
  cur_sigma2 = c()
  
  free_coef = c()
  free_cov = list()
  
  cov_list = list()
  
  split_labels = cut(x = Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(Data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  for(i in 1:n_dim){
    
    cur_regression = runUnboundedCubicSplineRegression_Path(Data,i,time_splits,split_labels, V_Smooth)
    
    if(is.null(cur_regression)){
      return(NULL)
    }
    
    cur_coef = cur_regression[[1]]
    
    free_coef = rbind(free_coef, cur_coef)
    cur_sigma2 = c(cur_sigma2, cur_regression[[2]])
    free_cov[[i]] = cur_regression[[3]]
    
  }
  
  trans_mat = constructUnboundedCubicSplineTransitionMatrix_Path(time_splits, V_Smooth)
  
  full_coef = t(trans_mat %*% t(free_coef))
  
  EstPath = lapply(unname(split(seq_len(ncol(full_coef)), ceiling(seq_len(ncol(full_coef)) / 4))), function(cols) full_coef[, cols, drop = FALSE])
  
  pos_splits = unname(table(split_labels_unique) >= 8)
  
  est_pos = getEstPathPosition(t = Data$t, EstPath = list(Path = EstPath, TimeSplits = time_splits))
  
  Error_mat = cbind(t = Data$t,(Data[,-1] - est_pos)^2)

  cur_full_cov_list = lapply(free_cov, FUN = function(x){trans_mat %*% x %*% t(trans_mat)})
  
  restricted_cov_perdim = lapply(cur_full_cov_list, FUN = function(cov){lapply(as.list(1:(length(time_splits) - 1)), FUN = function(x){cov[1:4 + 4*(x-1),1:4 + 4*(x-1)]})})
  cov_list = lapply(as.list(1:(length(time_splits)-1)), FUN = function(x){lapply(restricted_cov_perdim, `[[`, x)})
  sigma2_list = replicate(length(time_splits)-1, list(cur_sigma2))
  
  list(Path = EstPath, TimeSplits = time_splits, Error = Error_mat, n_splits = length(time_splits)-1, PosSplits = pos_splits, Cov = cov_list, Sigma2 = sigma2_list, SplitLabels = split_labels, FullCov = cur_full_cov_list, Free_Parameters = free_coef, Free_Parameter_Cov = free_cov)
  
}

#' getEstimatedPaths_FullRegression
#' 
#' @description
#' One of the functions to estimate a path. FullRegression method only fits all segments during each iteration. 
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#' @param Model_Selection_Type a string. The model criterion. 
#' Must pick from the following: BIC, BICVel, AIC, VarPos, VarVel, VarPosVel. 
#'
#' @returns an EstimatedPath object.
#' @export
#' 
#' @examples
#' getEstimatedPaths_FullRegression(ExData, 100, FALSE, "BIC")
#'
getEstimatedPaths_FullRegression = function(Data, nsplits, V_Smooth = F, Model_Selection_Type = "BIC"){
  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))
  
  time_splits = stats::fivenum(Data$t)
  
  all_splits_est_path = list()
  
  all_splits_est_path[[1]] = runEstPathStep_FullRegression(Data, time_splits, V_Smooth)
  
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  if(sum(all_splits_est_path[[1]]$PosSplits) == 0){
    all_splits_est_path
    
  } else{
    
    while(length(time_splits)-1 <= nsplits){
      pre_len = length(all_splits_est_path)
      worst_split = which.max(unname(by(rowSums(all_splits_est_path[[length(all_splits_est_path)]]$Error[,-1]), all_splits_est_path[[length(all_splits_est_path)]]$SplitLabels, sum))*all_splits_est_path[[length(all_splits_est_path)]]$PosSplits)
      time_splits = c(time_splits[1:worst_split], stats::median(Data$t[all_splits_est_path[[length(all_splits_est_path)]]$SplitLabels == worst_split]), time_splits[(worst_split+1):length(time_splits)])
      
      EstPath_Step = tryCatch(expr = runEstPathStep_FullRegression(Data, time_splits, V_Smooth), 
                              error = function(e){NULL})
      
      if(is.null(EstPath_Step)){
        
        all_splits_est_path[[length(all_splits_est_path)+1]] = all_splits_est_path[[length(all_splits_est_path)]]
        
        all_splits_est_path[[length(all_splits_est_path)]]$PosSplits[worst_split] = FALSE
        
      } else{
        
        all_splits_est_path[[length(all_splits_est_path)+1]] = EstPath_Step
        
      }
      
      post_len = length(all_splits_est_path)
      svMisc::progress(length(time_splits) - 2, max.value = nsplits, progress.bar = T, console = T)
      
      if(post_len == pre_len){
        break
      }
      
      if(sum(all_splits_est_path[[length(all_splits_est_path)]]$PosSplits) == 0){
        break
      }
      
    }
  }
  
  getBestSplitsPath(EstPath = all_splits_est_path, show = F, type = Model_Selection_Type, V_Smooth = V_Smooth)
  
}

