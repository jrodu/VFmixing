#' runEstPathStep_Projection
#' 
#' @description
#' At each iteration an unrestricted cubic function is fit to both halves of the splitting region. 
#' All unrestricted coefficients are projected into the restricted space.
#' 
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions. 
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param cur_unrestricted_coefs a list of matrices. Each element corresponds to the coefficient matrices for the unrestricted cubic functions fit on each region of the previous spline.
#' @param cur_unrestricted_cov a list of matrices. Each element corresponds to the covariance matrices for the unrestricted cubic coefficients fit on each region of the previous spline.
#' @param cur_integral_list a list of matrices. The previous integral list.
#' @param cur_sigma2_list a list of vectors. Each element corresponds to the sigma2 estimates of each dimension for the given region of the previous spline.
#' @param worst_split a numeric scalar. The region to be split.
#' @param V_Smooth a boolean. TRUE if the spline should have C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedPath object.
#' @export
#'
runEstPathStep_Projection = function(Data, time_splits, cur_unrestricted_coefs, cur_unrestricted_cov, cur_integral_list, cur_sigma2_list, worst_split, V_Smooth = T){
  
  n_dim = ncol(Data) - 1
  
  split_labels = cut(x = Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(Data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  first_coef = c()
  second_coef = c()
  
  first_Data = Data[split_labels == worst_split,]
  second_Data = Data[split_labels == worst_split+1,]
  
  tmp_unrestricted_cov = replicate(2, list())
  
  sigma2 = c()
  
  for(i in 1:n_dim){
    
    first_reg = runUnRestrictedCubicRegression_Path(cur_Data = first_Data, start_time = time_splits[worst_split], dim = i)
    second_reg = runUnRestrictedCubicRegression_Path(cur_Data = second_Data, start_time = time_splits[worst_split+1], dim = i)
    
    
    first_coef = rbind(first_coef, first_reg[[1]])
    second_coef = rbind(second_coef, second_reg[[1]])
    
    first_X = first_reg[[3]]
    second_X = second_reg[[3]]
    
    qr_first_XtX = qr(t(first_X) %*% first_X)
    
    if(qr_first_XtX$rank == 4){
      first_XtX_inv = backsolve(qr.R(qr_first_XtX), x = diag(4)) %*% t(qr.Q(qr_first_XtX))
    } else{
      first_XtX_inv = pracma::pinv(t(first_X) %*% first_X)
    }
    
    qr_second_XtX = qr(t(second_X) %*% second_X)
    
    if(qr_second_XtX$rank == 4){
      second_XtX_inv = backsolve(qr.R(qr_second_XtX), x = diag(4)) %*% t(qr.Q(qr_second_XtX))
    } else{
      second_XtX_inv = pracma::pinv(t(second_X) %*% second_X)
    }
    
    tmp_unrestricted_cov[[1]][[i]] = (first_reg[[2]] * first_XtX_inv)
    tmp_unrestricted_cov[[2]][[i]] = (second_reg[[2]] * second_XtX_inv)
    
    
  }
  
  cur_unrestricted_coefs[[worst_split]] = first_coef
  unrestricted_coefs_list = append(cur_unrestricted_coefs, list(second_coef), after = worst_split)
  unrestricted_coefs_mat = do.call("cbind", unrestricted_coefs_list)
  
  cur_unrestricted_cov[[worst_split]] = tmp_unrestricted_cov[[1]]
  unrestricted_cov = append(cur_unrestricted_cov, list(tmp_unrestricted_cov[[2]]), after = worst_split)
  
  
  integral_list = constructIntegralMatrixList_Path(cur_integral_list, time_splits, worst_split)
  
  Z = as.matrix(Matrix::bdiag(integral_list))
  
  V = constructUnboundedCubicSplineTransitionMatrix_Path(time_splits, V_Smooth = V_Smooth)
  
  VtZV = t(V) %*% Z %*% V
  
  qr_decomp_VtZV = qr(VtZV)
  
  if(qr_decomp_VtZV$rank == nrow(VtZV)){
    
    VtZV_inv = backsolve(qr.R(qr_decomp_VtZV), x = diag(ncol(V))) %*% t(qr.Q(qr_decomp_VtZV))
    
  } else{
    VtZV_inv = pracma::pinv(VtZV)
  }
  
  Free_Param_Matrix = VtZV_inv %*% t(V) %*% Z
  
  Hat_Matrix = (V %*% Free_Param_Matrix)
  
  free_param_coefs = t(Free_Param_Matrix %*% t(unrestricted_coefs_mat))
  
  restricted_coefs = t(Hat_Matrix %*% t(unrestricted_coefs_mat))
  
  EstPath = lapply(unname(split(seq_len(ncol(restricted_coefs)), ceiling(seq_len(ncol(restricted_coefs)) / 4))), function(cols) restricted_coefs[, cols, drop = FALSE])
  
  unrestricted_cov_perdim = lapply(as.list(1:n_dim), FUN = function(x){lapply(unrestricted_cov, `[[`, x)})
  
  full_restricted_cov = lapply(unrestricted_cov_perdim, FUN = function(cov){Hat_Matrix %*% as.matrix(Matrix::bdiag(cov)) %*% t(Hat_Matrix)})
  
  restricted_cov_perdim = lapply(full_restricted_cov, FUN = function(cov){lapply(as.list(1:(length(time_splits) - 1)), FUN = function(x){cov[1:4 + 4*(x-1),1:4 + 4*(x-1)]})})
  
  restricted_cov = lapply(as.list(1:(length(time_splits)-1)), FUN = function(x){lapply(restricted_cov_perdim, `[[`, x)})
  
  free_param_cov = lapply(unrestricted_cov_perdim, FUN = function(cov){Free_Param_Matrix %*% as.matrix(Matrix::bdiag(cov)) %*% t(Free_Param_Matrix)})
  
  est_pos = c()
  
  for(i in 1:(length(time_splits)-1)){
    est_pos = rbind(est_pos, evaluateCubic(Data$t[split_labels == i],Cubic = EstPath[[i]], startTime = time_splits[i]))
  }
  
  Error_mat = cbind(t = Data$t,(Data[,-1] - est_pos)^2)
  
  Sigma2_list = replicate(length(time_splits)-1, list(unname(colMeans(Error_mat[,-1]))))
  
  pos_splits = unname(table(split_labels_unique) >= 8)
  
  list(Path = EstPath, TimeSplits = time_splits, Error = Error_mat, n_splits = length(time_splits)-1, PosSplits = pos_splits, Cov = restricted_cov, FullCov = full_restricted_cov, Sigma2 = Sigma2_list, SplitLabels = split_labels, UnrestrictedCoef = unrestricted_coefs_list, UnrestrictedCov = unrestricted_cov, IntegralList = integral_list, Free_Parameters = free_param_coefs, Free_Parameter_Cov = free_param_cov)
  
}


#' getEstimatedPaths_Projection
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
#' getEstimatedPaths_Projection(ExData, 100, V_Smooth = FALSE, Model_Selection_Type = "BIC") 
#' 
getEstimatedPaths_Projection = function(Data, nsplits, V_Smooth = T, Model_Selection_Type = "VarPosVel"){
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))
  
  time_splits = c(min(Data$t),stats::median(Data$t),max(Data$t))
  
  all_splits_est_path = list()
  
  all_splits_est_path[[1]] = runEstPathStep_Projection(Data = Data, time_splits = time_splits, cur_unrestricted_coefs = list(), cur_unrestricted_cov = list(list()), cur_integral_list = list(), cur_sigma2_list = list(), worst_split = 1, V_Smooth = V_Smooth)
  
  svMisc::progress(1, max.value = nsplits, progress.bar = T, console = T)
  
  if(sum(all_splits_est_path[[1]]$PosSplits) == 0){
    all_splits_est_path
    
  } else{
    
    while(length(time_splits)-1 <= nsplits){
      
      worst_split = which.max(unname(by(rowSums(all_splits_est_path[[length(all_splits_est_path)]]$Error[,-1]), all_splits_est_path[[length(all_splits_est_path)]]$SplitLabels, sum))*all_splits_est_path[[length(all_splits_est_path)]]$PosSplits)
      time_splits = c(time_splits[1:worst_split], stats::median(Data$t[all_splits_est_path[[length(all_splits_est_path)]]$SplitLabels == worst_split]), time_splits[(worst_split+1):length(time_splits)])
      
      all_splits_est_path[[length(all_splits_est_path)+1]] = runEstPathStep_Projection(Data = Data, time_splits = time_splits, cur_unrestricted_coefs = all_splits_est_path[[length(all_splits_est_path)]]$UnrestrictedCoef, cur_unrestricted_cov = all_splits_est_path[[length(all_splits_est_path)]]$UnrestrictedCov, cur_integral_list = all_splits_est_path[[length(all_splits_est_path)]]$IntegralList, cur_sigma2_list = all_splits_est_path[[length(all_splits_est_path)]]$Sigma2, worst_split = worst_split, V_Smooth = V_Smooth)
      
      svMisc::progress(length(time_splits) - 2, max.value = nsplits, progress.bar = T, console = T)
      
      if(sum(all_splits_est_path[[length(all_splits_est_path)]]$PosSplits) == 0){
        break
      }
      
    }
  }
  
  getBestSplitsPath(EstPath = all_splits_est_path, show = F, type = Model_Selection_Type, V_Smooth = V_Smooth)
}



