#' getEstimatedTrajectories_ByFours
#' 
#' @description
#' By fours method for estimating the trajectory functions.
#' 
#'
#' @param EstPath an EstimatedPath object.
#' @param EstVel an EstimatedPath object. A cubic spline could be fit to estimates of velocity from the data. This is not required.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#' @param baseVectorFields a baseVectorField function.
#' @param data_int a numeric scalar. The time grid size used to constrict the trajectory data set.
#' @param Random_Path a boolean. TRUE if the path used to construct the trajectory data set should be drawn at random. 
#' If TRUE, the path must've been estimated using the FullRegression or Projection methods.
#' @param Path_V_Smooth a boolean. TRUE if the path cubic spline was estimated under C2 smoothness. FALSE if C1. Only needed if Random_Path = TRUE.
#'
#' @returns an EstimatedTrajectory object.
#' @export
#'
getEstimatedTrajectories_ByFours = function(EstPath, EstVel, nsplits, baseVectorFields, data_int = 0.001, Random_Path = F, Path_V_Smooth){
  
  n_dim = nrow(EstPath$Path[[1]])
  
  t_seq = seq(min(EstPath$TimeSplits), max(EstPath$TimeSplits), by = data_int)
  
  if(Random_Path){
    
    if(missing(EstVel)){
      rand_path = list(Path = generateRandomPath(EstPath$Free_Parameters, EstPath$Free_Parameter_Cov, EstPath$TimeSplits, Path_V_Smooth), TimeSplits = EstPath$TimeSplits)
      path_data = cbind(t_seq, getEstPathPosition(t_seq, rand_path), getEstPathVelocity(t_seq, rand_path))
      
    } else{
      rand_path = list(Path = generateRandomPath(EstPath$Free_Parameters, EstPath$Free_Parameter_Cov, EstPath$TimeSplits, Path_V_Smooth), TimeSplits = EstPath$TimeSplits)
      rand_vel = list(Path = generateRandomPath(EstVel$Free_Parameters, EstVel$Free_Parameter_Cov, EstVel$TimeSplits, Path_V_Smooth), TimeSplits = EstVel$TimeSplits)
      
      path_data = cbind(t_seq, getEstPathPosition(t_seq, rand_path), getEstPathPosition(t_seq, rand_vel))
    }
  } else{
    if(missing(EstVel)){
      path_data = cbind(t_seq, getEstPathPosition(t_seq, EstPath), getEstPathVelocity(t_seq, EstPath))
    } else{
      path_data = cbind(t_seq, getEstPathPosition(t_seq, EstPath), getEstPathPosition(t_seq, EstVel))
    }
  }
  
  colnames(path_data) = c('t', stringr::str_c("X", 1:n_dim), stringr::str_c("X", 1:n_dim, 'v'))
  
  base_vel = t(apply(path_data[,1:(n_dim+1)], MARGIN = 1, FUN = flatVF, baseVectorFields = baseVectorFields, simplify = T))
  
  not_0 = which(rowSums(base_vel == 0) == 0)
  
  n_models = ncol(base_vel)/n_dim
  
  colnames(base_vel) = c(stringr::str_c("X", rep(1:n_dim, n_models),'v', 'M',rep(1:n_models, each = n_dim)))
  
  
  path_data = cbind(path_data, base_vel)
  path_data = data.frame(path_data)[not_0,]
  
  all_splits_est_traj = list()
  
  # 0 Splits
  
  time_splits = unname(stats::quantile(path_data$t, probs = c(0,1)))
  split_labels = cut(x = path_data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(path_data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  est_traj = list()
  
  cur_traj1 = c()
  cur_Data = path_data
  
  cur_X = c()
  cur_Y = c()
  
  for(i in 1:n_dim){
    
    tmp_X = c()
    
    for(j in 1:n_models){
      
      cur_a1 = ((cur_Data$t - time_splits[1])^3)
      cur_b1 = ((cur_Data$t - time_splits[1])^2)
      cur_c1 = cur_Data$t - time_splits[1]
      cur_d1 = 1*(split_labels == 1)
      
      tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a1,cur_b1,cur_c1,cur_d1))
      
    }
    
    cur_X = rbind(cur_X, tmp_X)
    cur_Y = c(cur_Y, cur_Data[,1+n_dim+i])
    
  }
  
  cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
  
  for(i in 1:n_models){
    
    cur_traj1 = rbind(cur_traj1, unname(cur_coef[i,1:4]))
    
  }
  
  est_traj[[1]] = cur_traj1
  
  
  ### Calculating Initial Error Matrix
  
  cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]))
  
  est_vel = c()
  
  for(j in 1:n_dim){
    est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
  }
  
  Error_mat = cbind(t = cur_Data$t, (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2)
  
  pos_splits = unname(table(split_labels_unique) >= 8*n_models/n_dim)
  
  all_splits_est_traj[[1]] = list(Traj = est_traj, TimeSplits = time_splits, Error = Error_mat, n_splits = 0, PosSplits = pos_splits)
  
  ### 1 Split
  
  time_splits = unname(stats::quantile(path_data$t, probs = c(0,0.5,1)))
  split_labels = cut(x = path_data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(path_data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  est_traj = list()
  
  cur_traj1 = c()
  cur_traj2 = c()
  cur_Data = path_data
  
  cur_X = c()
  cur_Y = c()
  
  for(i in 1:n_dim){
    
    tmp_X = c()
    
    for(j in 1:n_models){
      
      cur_a1 = ((cur_Data$t - time_splits[1])^3)*(split_labels == 1) +
        (3*(time_splits[2]-time_splits[1])^2*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^3)*(split_labels %in% c(2:4))
      cur_b1 = ((cur_Data$t - time_splits[1])^2)*(split_labels == 1) +
        (2*(time_splits[2]-time_splits[1])*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^2)*(split_labels %in% c(2:4))
      cur_c1 = cur_Data$t - time_splits[1]
      cur_d1 = 1*(split_labels %in% c(1:4))
      cur_a2 = ((cur_Data$t - time_splits[2])^3)*(split_labels == 2)
      cur_b2 = ((cur_Data$t - time_splits[2])^2)*(split_labels == 2)
      
      tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a1,cur_b1,cur_c1,cur_d1,cur_a2,cur_b2))
      
    }
    
    cur_X = rbind(cur_X, tmp_X)
    cur_Y = c(cur_Y, cur_Data[,1+n_dim+i])
    
  }
  
  cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
  
  for(i in 1:n_models){
    
    cur_traj1 = rbind(cur_traj1, unname(cur_coef[i,1:4]))
    cur_traj2 = rbind(cur_traj2, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[2], cur_traj1[i,], startTime = time_splits[1]), evaluateCubic(time_splits[2], cur_traj1[i,], startTime = time_splits[1])))
    
  }
  
  est_traj[[1]] = cur_traj1
  est_traj[[2]] = cur_traj2
  
  
  ### Calculating Initial Error Matrix
  
  cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]),
                         evaluateCubic(path_data$t[split_labels == 2],Cubic = est_traj[[2]], startTime = time_splits[2]))
  
  est_vel = c()
  
  for(j in 1:n_dim){
    est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
  }
  
  Error_mat = cbind(t = cur_Data$t, (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2)
  
  pos_splits = unname(table(split_labels_unique) >= 8*n_models/n_dim)
  
  all_splits_est_traj[[2]] = list(Traj = est_traj, TimeSplits = time_splits, Error = Error_mat, n_splits = 1, PosSplits = pos_splits)
  
  # 2 Splits 
  
  time_splits = unname(stats::quantile(path_data$t, probs = c(0,1/3,2/3,1)))
  split_labels = cut(x = path_data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(path_data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  est_traj = list()
  
  cur_traj1 = c()
  cur_traj2 = c()
  cur_traj3 = c()
  cur_Data = path_data
  
  cur_X = c()
  cur_Y = c()
  
  for(i in 1:n_dim){
    
    tmp_X = c()
    
    for(j in 1:n_models){
      
      cur_a1 = ((cur_Data$t - time_splits[1])^3)*(split_labels == 1) +
        (3*(time_splits[2]-time_splits[1])^2*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^3)*(split_labels %in% c(2:4))
      cur_b1 = ((cur_Data$t - time_splits[1])^2)*(split_labels == 1) +
        (2*(time_splits[2]-time_splits[1])*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^2)*(split_labels %in% c(2:4))
      cur_c1 = cur_Data$t - time_splits[1]
      cur_d1 = 1*(split_labels %in% c(1:4))
      cur_a2 = ((cur_Data$t - time_splits[2])^3)*(split_labels == 2) +
        (3*(time_splits[3]-time_splits[2])^2*(cur_Data$t - time_splits[3]) + (time_splits[3] - time_splits[2])^3)*(split_labels %in% c(3:4))
      cur_b2 = ((cur_Data$t - time_splits[2])^2)*(split_labels == 2) +
        (2*(time_splits[3]-time_splits[2])*(cur_Data$t - time_splits[3]) + (time_splits[3] - time_splits[2])^2)*(split_labels %in% c(3:4))
      cur_a3 = ((cur_Data$t - time_splits[3])^3)*(split_labels == 3)
      cur_b3 = ((cur_Data$t - time_splits[3])^2)*(split_labels == 3)
      
      tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a1,cur_b1,cur_c1,cur_d1,cur_a2,cur_b2,cur_a3,cur_b3))
      
    }
    
    cur_X = rbind(cur_X, tmp_X)
    cur_Y = c(cur_Y, cur_Data[,1+n_dim+i])
    
  }
  
  cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
  
  for(i in 1:n_models){
    
    cur_traj1 = rbind(cur_traj1, unname(cur_coef[i,1:4]))
    cur_traj2 = rbind(cur_traj2, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[2], cur_traj1[i,], startTime = time_splits[1]), evaluateCubic(time_splits[2], cur_traj1[i,], startTime = time_splits[1])))
    cur_traj3 = rbind(cur_traj3, c(unname(cur_coef[i,7:8]), evaluateCubicVelocity(time_splits[3], cur_traj2[i,], startTime = time_splits[2]), evaluateCubic(time_splits[3], cur_traj2[i,], startTime = time_splits[2])))
    
  }
  
  est_traj[[1]] = cur_traj1
  est_traj[[2]] = cur_traj2
  est_traj[[3]] = cur_traj3
  
  
  ### Calculating Initial Error Matrix
  
  cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]),
                         evaluateCubic(path_data$t[split_labels == 2],Cubic = est_traj[[2]], startTime = time_splits[2]),
                         evaluateCubic(path_data$t[split_labels == 3],Cubic = est_traj[[3]], startTime = time_splits[3]))
  
  est_vel = c()
  
  for(j in 1:n_dim){
    est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
  }
  
  Error_mat = cbind(t = cur_Data$t, (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2)
  
  pos_splits = unname(table(split_labels_unique) >= 8*n_models/n_dim)
  
  all_splits_est_traj[[3]] = list(Traj = est_traj, TimeSplits = time_splits, Error = Error_mat, n_splits = 2, PosSplits = pos_splits)
  
  ## 4 Splits
  
  time_splits = stats::fivenum(path_data$t)
  split_labels = cut(x = path_data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(path_data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  est_traj = list()
  
  cur_traj1 = c()
  cur_traj2 = c()
  cur_traj3 = c()
  cur_traj4 = c()
  cur_Data = path_data
  
  cur_X = c()
  cur_Y = c()
  
  for(i in 1:n_dim){
    
    tmp_X = c()
    
    for(j in 1:n_models){
      
      cur_a1 = ((cur_Data$t - time_splits[1])^3)*(split_labels == 1) +
        (3*(time_splits[2]-time_splits[1])^2*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^3)*(split_labels %in% c(2:4))
      cur_b1 = ((cur_Data$t - time_splits[1])^2)*(split_labels == 1) +
        (2*(time_splits[2]-time_splits[1])*(cur_Data$t - time_splits[2]) + (time_splits[2] - time_splits[1])^2)*(split_labels %in% c(2:4))
      cur_c1 = cur_Data$t - time_splits[1]
      cur_d1 = 1*(split_labels %in% c(1:4))
      cur_a2 = ((cur_Data$t - time_splits[2])^3)*(split_labels == 2) +
        (3*(time_splits[3]-time_splits[2])^2*(cur_Data$t - time_splits[3]) + (time_splits[3] - time_splits[2])^3)*(split_labels %in% c(3:4))
      cur_b2 = ((cur_Data$t - time_splits[2])^2)*(split_labels == 2) +
        (2*(time_splits[3]-time_splits[2])*(cur_Data$t - time_splits[3]) + (time_splits[3] - time_splits[2])^2)*(split_labels %in% c(3:4))
      cur_a3 = ((cur_Data$t - time_splits[3])^3)*(split_labels == 3) +
        (3*(time_splits[4]-time_splits[3])^2*(cur_Data$t - time_splits[4]) + (time_splits[4] - time_splits[3])^3)*(split_labels == 4)
      cur_b3 = ((cur_Data$t - time_splits[3])^2)*(split_labels == 3) +
        (2*(time_splits[4]-time_splits[3])*(cur_Data$t - time_splits[4]) + (time_splits[4] - time_splits[3])^2)*(split_labels == 4)
      cur_a4 = ((cur_Data$t - time_splits[4])^3)*(split_labels == 4)
      cur_b4 = ((cur_Data$t - time_splits[4])^2)*(split_labels == 4)
      
      tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a1,cur_b1,cur_c1,cur_d1,cur_a2,cur_b2,cur_a3,cur_b3,cur_a4,cur_b4))
      
    }
    
    cur_X = rbind(cur_X, tmp_X)
    cur_Y = c(cur_Y, cur_Data[,1+n_dim+i])
    
  }
  
  cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
  
  for(i in 1:n_models){
    
    cur_traj1 = rbind(cur_traj1, unname(cur_coef[i,1:4]))
    cur_traj2 = rbind(cur_traj2, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[2], cur_traj1[i,], startTime = time_splits[1]), evaluateCubic(time_splits[2], cur_traj1[i,], startTime = time_splits[1])))
    cur_traj3 = rbind(cur_traj3, c(unname(cur_coef[i,7:8]), evaluateCubicVelocity(time_splits[3], cur_traj2[i,], startTime = time_splits[2]), evaluateCubic(time_splits[3], cur_traj2[i,], startTime = time_splits[2])))
    cur_traj4 = rbind(cur_traj4, c(unname(cur_coef[i,9:10]), evaluateCubicVelocity(time_splits[4], cur_traj3[i,], startTime = time_splits[3]), evaluateCubic(time_splits[4], cur_traj3[i,], startTime = time_splits[3])))
    
  }
  
  est_traj[[1]] = cur_traj1
  est_traj[[2]] = cur_traj2
  est_traj[[3]] = cur_traj3
  est_traj[[4]] = cur_traj4
  
  ### Calculating Initial Error Matrix
  
  cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]),
                         evaluateCubic(path_data$t[split_labels == 2],Cubic = est_traj[[2]], startTime = time_splits[2]),
                         evaluateCubic(path_data$t[split_labels == 3],Cubic = est_traj[[3]], startTime = time_splits[3]),
                         evaluateCubic(path_data$t[split_labels == 4],Cubic = est_traj[[4]], startTime = time_splits[4]))
  
  est_vel = c()
  
  for(j in 1:n_dim){
    est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
  }
  
  Error_mat = cbind(t = cur_Data$t, (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2)
  
  pos_splits = unname(table(split_labels_unique) >= 8*n_models/n_dim)
  
  all_splits_est_traj[[4]] = list(Traj = est_traj, TimeSplits = time_splits, Error = Error_mat, n_splits = 3, PosSplits = pos_splits)
  
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  while(length(time_splits)-1 <= nsplits){
    
    split_labels_path = cut(x = EstPath$Error$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
    mean_path_error = unname(by(rowSums(EstPath$Error[,-1]), split_labels_path, stats::median))
    mean_path_error[is.na(mean_path_error)] = stats::median(mean_path_error, na.rm = T)
    
    worst_split = which.max((log(unname(by(rowSums(all_splits_est_traj[[length(all_splits_est_traj)]]$Error[,-1]), split_labels, sum)+1)) - mean_path_error + max(mean_path_error))*pos_splits)
    
    time_splits = c(time_splits[1:worst_split], stats::median(path_data$t[split_labels == worst_split]), time_splits[(worst_split+1):length(time_splits)])
    split_labels = cut(x = path_data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
    split_labels_unique = cut(x = unique(path_data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
    
    if(worst_split == 1){
      
      cur_Data = path_data[split_labels %in% c(1:3),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_traj_first = c()
      cur_traj_second = c()
      cur_traj_next = c()
      cur_X = c()
      cur_Y = c()
      
      for(i in 1:n_dim){
        
        tmp_X = c()
        tmp_C = c()
        
        for(j in 1:n_models){
          
          end_pos = est_traj[[3]][j,4]
          end_vel = est_traj[[3]][j,3]
          
          for(k in 1:3){
            
            cur_a = (cur_Data$t - time_splits[k])^3*(cur_split_labels == k) +
              (2*(time_splits[k+1]-time_splits[k])^3-3*(time_splits[k+1]-time_splits[k])^2*(cur_Data$t - time_splits[k]))*(cur_split_labels %in% c(1:k))
            
            cur_b = (cur_Data$t - time_splits[k])^2*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^2-2*(time_splits[k+1]-time_splits[k])*(cur_Data$t - time_splits[k]))*(cur_split_labels %in% c(1:k))
            
            tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a, cur_b))
            
          }
          
          tmp_C = cbind(tmp_C, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*(end_pos + end_vel*(cur_Data$t - time_splits[4])))
          
        }
        
        cur_X = rbind(cur_X, tmp_X)
        cur_Y = c(cur_Y, cur_Data[,1+n_dim+i] - rowSums(tmp_C))
        
      }
      
      cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
      
      for(i in 1:n_models){
        
        end_pos = est_traj[[3]][i,4]
        end_vel = est_traj[[3]][i,3]
        
        c_next = end_vel - 3*cur_coef[i,5]*(time_splits[4]-time_splits[3])^2 - 2*cur_coef[i,6]*(time_splits[4]-time_splits[3])
        d_next = end_pos - cur_coef[i,5]*(time_splits[4]-time_splits[3])^3 - cur_coef[i,6]*(time_splits[4]-time_splits[3])^2 - c_next*(time_splits[4]-time_splits[3])
        cur_traj_next = rbind(cur_traj_next, c(unname(cur_coef[i,5:6]), c_next, d_next))
        
        c_second = cur_traj_next[i,3] - 3*cur_coef[i,3]*(time_splits[3]-time_splits[2])^2 - 2*cur_coef[i,4]*(time_splits[3]-time_splits[2])
        d_second = cur_traj_next[i,4] - cur_coef[i,3]*(time_splits[3]-time_splits[2])^3 - cur_coef[i,4]*(time_splits[3]-time_splits[2])^2 - c_second*(time_splits[3]-time_splits[2])
        cur_traj_second = rbind(cur_traj_second, c(unname(cur_coef[i,3:4]), c_second, d_second))
        
        c_first = cur_traj_second[i,3] - 3*cur_coef[i,1]*(time_splits[2]-time_splits[1])^2 - 2*cur_coef[i,2]*(time_splits[2]-time_splits[1])
        d_first = cur_traj_second[i,4] - cur_coef[i,1]*(time_splits[2]-time_splits[1])^3 - cur_coef[i,2]*(time_splits[2]-time_splits[1])^2 - c_first*(time_splits[2]-time_splits[1])
        cur_traj_first = rbind(cur_traj_first, c(unname(cur_coef[i,1:2]), c_first, d_first))
        
      }
      
      est_traj[[1]] = cur_traj_first
      est_traj[[2]] = cur_traj_next
      est_traj = append(est_traj, list(cur_traj_second), after = 1)
      
      cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]),
                             evaluateCubic(path_data$t[split_labels == 2],Cubic = est_traj[[2]], startTime = time_splits[2]),
                             evaluateCubic(path_data$t[split_labels == 3],Cubic = est_traj[[3]], startTime = time_splits[3]))
      
      est_vel = c()
      
      for(j in 1:n_dim){
        est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
      }
      
      cur_Error_mat = (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2
      
      Error_mat[split_labels %in% c(1:3),-1] = cur_Error_mat
      
    } else if(worst_split == 2){
      
      cur_Data = path_data[split_labels %in% c(1:4),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_traj_prev = c()
      cur_traj_first = c()
      cur_traj_second = c()
      cur_traj_next = c()
      cur_X = c()
      cur_Y = c()
      
      for(i in 1:n_dim){
        
        tmp_X = c()
        tmp_C = c()
        
        for(j in 1:n_models){
          
          end_pos = est_traj[[4]][j,4]
          end_vel = est_traj[[4]][j,3]
          
          for(k in 1:4){
            
            cur_a = (cur_Data$t - time_splits[k])^3*(cur_split_labels == k) +
              (2*(time_splits[k+1]-time_splits[k])^3-3*(time_splits[k+1]-time_splits[k])^2*(cur_Data$t - time_splits[k]))*(cur_split_labels %in% c(1:k))
            
            cur_b = (cur_Data$t - time_splits[k])^2*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^2-2*(time_splits[k+1]-time_splits[k])*(cur_Data$t - time_splits[k]))*(cur_split_labels %in% c(1:k))
            
            tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a, cur_b))
            
          }
          
          tmp_C = cbind(tmp_C, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*(end_pos + end_vel*(cur_Data$t - time_splits[5])))
          
        }
        
        cur_X = rbind(cur_X, tmp_X)
        cur_Y = c(cur_Y, cur_Data[,1+n_dim+i] - rowSums(tmp_C))
        
      }
      
      cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
      
      for(i in 1:n_models){
        
        end_pos = est_traj[[4]][i,4]
        end_vel = est_traj[[4]][i,3]
        
        c_next = end_vel - 3*cur_coef[i,7]*(time_splits[5]-time_splits[4])^2 - 2*cur_coef[i,8]*(time_splits[5]-time_splits[4])
        d_next = end_pos - cur_coef[i,7]*(time_splits[5]-time_splits[4])^3 - cur_coef[i,8]*(time_splits[5]-time_splits[4])^2 - c_next*(time_splits[5]-time_splits[4])
        cur_traj_next = rbind(cur_traj_next, c(unname(cur_coef[i,7:8]), c_next, d_next))
        
        c_second = cur_traj_next[i,3] - 3*cur_coef[i,5]*(time_splits[4]-time_splits[3])^2 - 2*cur_coef[i,6]*(time_splits[4]-time_splits[3])
        d_second = cur_traj_next[i,4] - cur_coef[i,5]*(time_splits[4]-time_splits[3])^3 - cur_coef[i,6]*(time_splits[4]-time_splits[3])^2 - c_second*(time_splits[4]-time_splits[3])
        cur_traj_second = rbind(cur_traj_second, c(unname(cur_coef[i,5:6]), c_second, d_second))
        
        c_first = cur_traj_second[i,3] - 3*cur_coef[i,3]*(time_splits[3]-time_splits[2])^2 - 2*cur_coef[i,4]*(time_splits[3]-time_splits[2])
        d_first = cur_traj_second[i,4] - cur_coef[i,3]*(time_splits[3]-time_splits[2])^3 - cur_coef[i,4]*(time_splits[3]-time_splits[2])^2 - c_first*(time_splits[3]-time_splits[2])
        cur_traj_first = rbind(cur_traj_first, c(unname(cur_coef[i,3:4]), c_first, d_first))
        
        c_prev = cur_traj_first[i,3] - 3*cur_coef[i,1]*(time_splits[2]-time_splits[1])^2 - 2*cur_coef[i,2]*(time_splits[2]-time_splits[1])
        d_prev = cur_traj_first[i,4] - cur_coef[i,1]*(time_splits[2]-time_splits[1])^3 - cur_coef[i,2]*(time_splits[2]-time_splits[1])^2 - c_prev*(time_splits[2]-time_splits[1])
        cur_traj_prev = rbind(cur_traj_prev, c(unname(cur_coef[i,1:2]), c_prev, d_prev))
        
      }
      
      est_traj[[1]] = cur_traj_prev
      est_traj[[2]] = cur_traj_first
      est_traj[[3]] = cur_traj_next
      est_traj = append(est_traj, list(cur_traj_second), after = 2)
      
      cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == 1],Cubic = est_traj[[1]], startTime = time_splits[1]),
                             evaluateCubic(path_data$t[split_labels == 2],Cubic = est_traj[[2]], startTime = time_splits[2]),
                             evaluateCubic(path_data$t[split_labels == 3],Cubic = est_traj[[3]], startTime = time_splits[3]),
                             evaluateCubic(path_data$t[split_labels == 4],Cubic = est_traj[[4]], startTime = time_splits[4]))
      
      est_vel = c()
      
      for(j in 1:n_dim){
        est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
      }
      
      cur_Error_mat = (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2
      Error_mat[split_labels %in% c(1:4),-1] = cur_Error_mat
      
    } else if(worst_split == length(time_splits) - 3){
      
      cur_Data = path_data[split_labels %in% c(-1:2 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_traj_prev = c()
      cur_traj_first = c()
      cur_traj_second = c()
      cur_traj_next = c()
      cur_X = c()
      cur_Y = c()
      
      for(i in 1:n_dim){
        
        tmp_X = c()
        tmp_C = c()
        
        for(j in 1:n_models){
          
          init_pos = est_traj[[worst_split - 1]][j,4]
          init_vel = est_traj[[worst_split - 1]][j,3]
          
          for(k in -1:2 + worst_split){
            
            cur_a = (cur_Data$t - time_splits[k])^3*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^3+3*(time_splits[k+1]-time_splits[k])^2*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(length(time_splits)-1)))*(k != (length(time_splits)-1))
            
            cur_b = (cur_Data$t - time_splits[k])^2*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^2+2*(time_splits[k+1]-time_splits[k])*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(length(time_splits)-1)))*(k != (length(time_splits)-1))
            
            tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a, cur_b))
            
          }
          
          tmp_C = cbind(tmp_C, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*(init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1])))
          
        }
        
        cur_X = rbind(cur_X, tmp_X)
        cur_Y = c(cur_Y, cur_Data[,1+n_dim+i] - rowSums(tmp_C))
        
      }
      
      cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
      
      for(i in 1:n_models){
        
        init_pos = est_traj[[worst_split - 1]][i,4]
        init_vel = est_traj[[worst_split - 1]][i,3]
        
        cur_traj_prev = rbind(cur_traj_prev, c(unname(cur_coef[i,1:2]), init_vel, init_pos))
        cur_traj_first = rbind(cur_traj_first, c(unname(cur_coef[i,3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1])))
        cur_traj_second = rbind(cur_traj_second, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split])))
        cur_traj_next = rbind(cur_traj_next, c(unname(cur_coef[i,7:8]), evaluateCubicVelocity(time_splits[worst_split + 2], cur_traj_second[i,], time_splits[worst_split + 1]), evaluateCubic(time_splits[worst_split + 2], cur_traj_second[i,], time_splits[worst_split + 1])))
        
      }
      
      est_traj[[worst_split - 1]] = cur_traj_prev
      est_traj[[worst_split]] = cur_traj_first
      est_traj[[worst_split + 1]] = cur_traj_next
      est_traj = append(est_traj, list(cur_traj_second), after = worst_split)
      
      cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == worst_split - 1],Cubic = est_traj[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                             evaluateCubic(path_data$t[split_labels == worst_split],Cubic = est_traj[[worst_split]], startTime = time_splits[worst_split]),
                             evaluateCubic(path_data$t[split_labels == worst_split + 1],Cubic = est_traj[[worst_split + 1]], startTime = time_splits[worst_split + 1]),
                             evaluateCubic(path_data$t[split_labels == worst_split + 2],Cubic = est_traj[[worst_split + 2]], startTime = time_splits[worst_split + 2]))
      
      est_vel = c()
      
      for(j in 1:n_dim){
        est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
      }
      
      cur_Error_mat = (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2
      
      Error_mat[split_labels %in% c(-1:2 + worst_split),-1] = cur_Error_mat
      
    } else if(worst_split == length(time_splits) - 2){
      
      cur_Data = path_data[split_labels %in% c(-1:1 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_traj_prev = c()
      cur_traj_first = c()
      cur_traj_second = c()
      cur_X = c()
      cur_Y = c()
      
      for(i in 1:n_dim){
        
        tmp_X = c()
        tmp_C = c()
        
        for(j in 1:n_models){
          
          init_pos = est_traj[[worst_split - 1]][j,4]
          init_vel = est_traj[[worst_split - 1]][j,3]
          
          for(k in -1:1 + worst_split){
            
            cur_a = (cur_Data$t - time_splits[k])^3*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^3+3*(time_splits[k+1]-time_splits[k])^2*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(length(time_splits)-1)))*(k != (length(time_splits)-1))
            
            cur_b = (cur_Data$t - time_splits[k])^2*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^2+2*(time_splits[k+1]-time_splits[k])*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(length(time_splits)-1)))*(k != (length(time_splits)-1))
            
            tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a, cur_b))
            
          }
          
          tmp_C = cbind(tmp_C, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*(init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1])))
          
        }
        
        cur_X = rbind(cur_X, tmp_X)
        cur_Y = c(cur_Y, cur_Data[,1+n_dim+i] - rowSums(tmp_C))
        
      }
      
      cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
      
      for(i in 1:n_models){
        
        init_pos = est_traj[[worst_split - 1]][i,4]
        init_vel = est_traj[[worst_split - 1]][i,3]
        
        cur_traj_prev = rbind(cur_traj_prev, c(unname(cur_coef[i,1:2]), init_vel, init_pos))
        cur_traj_first = rbind(cur_traj_first, c(unname(cur_coef[i,3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1])))
        cur_traj_second = rbind(cur_traj_second, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split])))
        
      }
      
      est_traj[[worst_split - 1]] = cur_traj_prev
      est_traj[[worst_split]] = cur_traj_first
      est_traj = append(est_traj, list(cur_traj_second), after = worst_split)
      
      cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == worst_split - 1],Cubic = est_traj[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                             evaluateCubic(path_data$t[split_labels == worst_split],Cubic = est_traj[[worst_split]], startTime = time_splits[worst_split]),
                             evaluateCubic(path_data$t[split_labels == worst_split + 1],Cubic = est_traj[[worst_split + 1]], startTime = time_splits[worst_split + 1]))
      
      est_vel = c()
      
      for(j in 1:n_dim){
        est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
      }
      
      cur_Error_mat = (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2
      
      Error_mat[split_labels %in% c(-1:1 + worst_split),-1] = cur_Error_mat
      
      
    } else{
      
      cur_Data = path_data[split_labels %in% c(-1:2 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_traj_prev = c()
      cur_traj_first = c()
      cur_traj_second = c()
      cur_traj_next = c()
      cur_X = c()
      cur_Y = c()
      
      for(i in 1:n_dim){
        
        tmp_X = c()
        tmp_C = c()
        
        for(j in 1:n_models){
          
          init_pos = est_traj[[worst_split - 1]][j,4]
          init_vel = est_traj[[worst_split - 1]][j,3]
          end_pos = est_traj[[worst_split + 2]][j,4]
          end_vel = est_traj[[worst_split + 2]][j,3]
          
          for(k in -1:1 + worst_split){
            
            cur_a = (cur_Data$t - time_splits[k])^3*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^3+3*(time_splits[k+1]-time_splits[k])^2*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(worst_split + 2))) +
              (((2*(time_splits[k+1] - time_splits[k])^3 + 6*(time_splits[k+1] - time_splits[k])^2*(time_splits[worst_split+2] - time_splits[k+1]) + 3*(time_splits[k+1] - time_splits[k])^2*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^3))*(cur_Data$t - time_splits[worst_split+2])^3 -
                 ((3*(time_splits[k+1] - time_splits[k])^3 + 9*(time_splits[k+1] - time_splits[k])^2*(time_splits[worst_split+2] - time_splits[k+1]) + 6*(time_splits[k+1] - time_splits[k])^2*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^2))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2)
            
            cur_b = (cur_Data$t - time_splits[k])^2*(cur_split_labels == k) +
              ((time_splits[k+1]-time_splits[k])^2+2*(time_splits[k+1]-time_splits[k])*(cur_Data$t - time_splits[k+1]))*(cur_split_labels %in% c((k+1):(worst_split + 2))) +
              (((2*(time_splits[k+1] - time_splits[k])^2 + 4*(time_splits[k+1] - time_splits[k])*(time_splits[worst_split+2] - time_splits[k+1]) + 2*(time_splits[k+1] - time_splits[k])*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^3))*(cur_Data$t - time_splits[worst_split+2])^3 -
                 ((3*(time_splits[k+1] - time_splits[k])^2 + 6*(time_splits[k+1] - time_splits[k])*(time_splits[worst_split+2] - time_splits[k+1]) + 4*(time_splits[k+1] - time_splits[k])*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^2))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2)
            
            tmp_X = cbind(tmp_X, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*cbind(cur_a, cur_b))
            
          }
          
          tmp_C = cbind(tmp_C, cur_Data[,1 + 2*n_dim + n_dim*(j-1) + i]*(init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1]) +
                                                                           (((end_vel + init_vel)/(time_splits[worst_split+3] - time_splits[worst_split+2])^2 + (2*init_vel*(time_splits[worst_split+2] - time_splits[worst_split-1]) - 2*(end_pos-init_pos))/(time_splits[worst_split+3] - time_splits[worst_split+2])^3)*(cur_Data$t - time_splits[worst_split+2])^3 +
                                                                              ((3*(end_pos - init_pos) - 3*init_vel*(time_splits[worst_split+2] - time_splits[worst_split-1]))/(time_splits[worst_split+3] - time_splits[worst_split+2])^2 - (end_vel + 2*init_vel)/(time_splits[worst_split+3] - time_splits[worst_split+2]))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2)))
        }
        
        cur_X = rbind(cur_X, tmp_X)
        cur_Y = c(cur_Y, cur_Data[,1+n_dim+i] - rowSums(tmp_C))
        
      }
      
      cur_coef = matrix(unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients), nrow = n_models, byrow = T)
      
      for(i in 1:n_models){
        
        init_pos = est_traj[[worst_split - 1]][i,4]
        init_vel = est_traj[[worst_split - 1]][i,3]
        end_pos = est_traj[[worst_split + 2]][i,4]
        end_vel = est_traj[[worst_split + 2]][i,3]
        
        cur_traj_prev = rbind(cur_traj_prev, c(unname(cur_coef[i,1:2]), init_vel, init_pos))
        cur_traj_first = rbind(cur_traj_first, c(unname(cur_coef[i,3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_traj_prev[i,], time_splits[worst_split - 1])))
        cur_traj_second = rbind(cur_traj_second, c(unname(cur_coef[i,5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_traj_first[i,], time_splits[worst_split])))
        
        init_vel_next = evaluateCubicVelocity(time_splits[worst_split + 2], cur_traj_second[i,], time_splits[worst_split+1])
        init_pos_next = evaluateCubic(time_splits[worst_split + 2], cur_traj_second[i,], time_splits[worst_split+1])
        diffTime_next = time_splits[worst_split+3] - time_splits[worst_split+2]
        
        cur_traj_next = rbind(cur_traj_next, c((init_vel_next + end_vel)/(diffTime_next^2) - 2*(end_pos - init_pos_next)/(diffTime_next^3),
                                               (end_vel - init_vel_next)/(2*diffTime_next) - 3*(init_vel_next + end_vel)/(2*diffTime_next) + 3*(end_pos - init_pos_next)/(diffTime_next^2),
                                               init_vel_next, init_pos_next))
        
      }
      
      est_traj[[worst_split - 1]] = cur_traj_prev
      est_traj[[worst_split]] = cur_traj_first
      est_traj[[worst_split + 1]] = cur_traj_next
      est_traj = append(est_traj, list(cur_traj_second), after = worst_split)
      
      
      cur_traj_value = rbind(evaluateCubic(path_data$t[split_labels == worst_split - 1],Cubic = est_traj[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                             evaluateCubic(path_data$t[split_labels == worst_split],Cubic = est_traj[[worst_split]], startTime = time_splits[worst_split]),
                             evaluateCubic(path_data$t[split_labels == worst_split + 1],Cubic = est_traj[[worst_split + 1]], startTime = time_splits[worst_split + 1]),
                             evaluateCubic(path_data$t[split_labels == worst_split + 2],Cubic = est_traj[[worst_split + 2]], startTime = time_splits[worst_split + 2]))
      
      
      est_vel = c()
      
      for(j in 1:n_dim){
        est_vel = cbind(est_vel, rowSums(cur_traj_value * as.matrix(cur_Data[,n_dim*(0:(n_models-1)) + 1 + 2*n_dim + j], ncol = n_models)))
      }
      
      cur_Error_mat = (cur_Data[,1:n_dim+1+n_dim] - est_vel)^2
      
      Error_mat[split_labels %in% c(-1:2 + worst_split),-1] = cur_Error_mat
      
    }
    
    pos_splits = unname(table(split_labels_unique) >= 8*n_models/n_dim)
    
    all_splits_est_traj[[length(all_splits_est_traj)+1]] = list(Traj = est_traj, TimeSplits = time_splits, Error = Error_mat, n_splits = length(est_traj)-1, PosSplits = pos_splits)
    
    svMisc::progress(length(est_traj)-1, max.value = nsplits, progress.bar = T, console = T)
    
    if(sum(pos_splits) == 0){
      break
    }
    
  }
  
  
  getBestSplitsTraj(all_splits_est_traj, EstPath = EstPath, show = F, type = "BIC", V_Smooth = F)
  
}
