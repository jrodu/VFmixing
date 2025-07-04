#' Get Estimated Path using the By Fours Method
#' 
#' @description
#' One of the functions to estimate a path. ByFours method only fits the previous segment, the two splitting segments, and the next segment during each iteration. 
#' 
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#'
#' @returns an EstimatedPath object. 
#' @export
#'
getEstimatedPaths_ByFours = function(Data, nsplits){
  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))
  
  time_splits = stats::fivenum(Data$t)
  
  all_splits_est_path = list()
  
  cur_path1 = c()
  cur_path2 = c()
  cur_path3 = c()
  cur_path4 = c()
  cur_Data = Data
  
  split_labels = cut(x = Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  split_labels_unique = cut(x = unique(Data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  
  runRegressionInit = function(dim,start_split, end_split,split_labels){
    
    a1 = ((cur_Data$t - time_splits[start_split])^3)*(split_labels == start_split) +
      (3*(time_splits[start_split+1]-time_splits[start_split])^2*(cur_Data$t - time_splits[start_split+1]) + (time_splits[start_split+1] - time_splits[start_split])^3)*(split_labels %in% c((start_split+1):end_split))
    b1 = ((cur_Data$t - time_splits[start_split])^2)*(split_labels == start_split) +
      (2*(time_splits[start_split+1]-time_splits[start_split])*(cur_Data$t - time_splits[start_split+1]) + (time_splits[start_split+1] - time_splits[start_split])^2)*(split_labels %in% c((start_split+1):end_split))
    c1 = cur_Data$t - time_splits[start_split]
    d1 = 1*(split_labels %in% c(start_split:end_split))
    
    cur_X = cbind(a1,b1,c1,d1)
    
    for(j in (start_split+1):end_split){
      cur_a = ((cur_Data$t - time_splits[j])^3)*(split_labels == j) +
        (3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j+1]) + (time_splits[j+1] - time_splits[j])^3)*(split_labels %in% c((j+1):end_split))*(j!=end_split)
      cur_b = ((cur_Data$t - time_splits[j])^2)*(split_labels == j) +
        (2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j+1]) + (time_splits[j+1] - time_splits[j])^2)*(split_labels %in% c((j+1):end_split))*(j!=end_split)
      cur_X = cbind(cur_X, cur_a, cur_b)
    }
    cur_Y = cur_Data[,dim+1]
    
    cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
    
    cur_coef
    
  }
  
  for(i in 1:n_dim){
    
    cur_coef = runRegressionInit(i,1,4,split_labels)
    
    cur_path1 = rbind(cur_path1, unname(cur_coef[1:4]))
    cur_path2 = rbind(cur_path2, c(unname(cur_coef[5:6]), evaluateCubicVelocity(time_splits[2], cur_path1[i,], startTime = time_splits[1]), evaluateCubic(time_splits[2], cur_path1[i,], startTime = time_splits[1])))
    cur_path3 = rbind(cur_path3, c(unname(cur_coef[7:8]), evaluateCubicVelocity(time_splits[3], cur_path2[i,], startTime = time_splits[2]), evaluateCubic(time_splits[3], cur_path2[i,], startTime = time_splits[2])))
    cur_path4 = rbind(cur_path4, c(unname(cur_coef[9:10]), evaluateCubicVelocity(time_splits[4], cur_path3[i,], startTime = time_splits[3]), evaluateCubic(time_splits[4], cur_path3[i,], startTime = time_splits[3])))
    
  }
  
  est_path[[1]] = cur_path1
  est_path[[2]] = cur_path2
  est_path[[3]] = cur_path3
  est_path[[4]] = cur_path4
  
  pos_splits = unname(table(split_labels_unique) >= 8)
  
  #Calculating Initial Error Matrix
  
  est_pos = rbind(evaluateCubic(Data$t[split_labels == 1],Cubic = est_path[[1]], startTime = time_splits[1]),
                  evaluateCubic(Data$t[split_labels == 2],Cubic = est_path[[2]], startTime = time_splits[2]),
                  evaluateCubic(Data$t[split_labels == 3],Cubic = est_path[[3]], startTime = time_splits[3]),
                  evaluateCubic(Data$t[split_labels == 4],Cubic = est_path[[4]], startTime = time_splits[4]))
  
  Error_mat = cbind(t = cur_Data$t,(cur_Data[,-1] - est_pos)^2)
  
  all_splits_est_path[[1]] = list(Path = est_path, TimeSplits = time_splits, Error = Error_mat, n_splits = 3, PosSplits = pos_splits)
  
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  while(length(time_splits)-1 <= nsplits){
    
    worst_split = which.max(unname(by(rowSums(all_splits_est_path[[length(all_splits_est_path)]]$Error[,-1]), split_labels, sum))*pos_splits)
    time_splits = c(time_splits[1:worst_split], stats::median(Data$t[split_labels == worst_split]), time_splits[(worst_split+1):length(time_splits)])
    split_labels = cut(x = Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
    split_labels_unique = cut(x = unique(Data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
    
    if(worst_split == 1){
      
      cur_Data = Data[split_labels %in% c(1:3),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_path_first = c()
      cur_path_second = c()
      cur_path_next = c()
      
      runRegressionWS1 = function(end_cond, dim){
        
        end_vel = end_cond[1]
        end_pos = end_cond[2]
        
        cur_X = c()
        
        for(j in 1:3){
          cur_a = (cur_Data$t - time_splits[j])^3*(cur_split_labels == j) +
            (2*(time_splits[j+1]-time_splits[j])^3-3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j]))*(cur_split_labels %in% c(1:j))
          
          cur_b = (cur_Data$t - time_splits[j])^2*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^2-2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j]))*(cur_split_labels %in% c(1:j))
          
          cur_X = cbind(cur_X, cur_a, cur_b)
          
        }
        
        cur_Y = cur_Data[,dim+1] - (end_pos + end_vel*(cur_Data$t - time_splits[4]))
        
        cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
        
        cur_coef
        
      }
      
      for(i in 1:n_dim){
        
        end_pos = est_path[[3]][i,4]
        end_vel = est_path[[3]][i,3]
        
        cur_coef = runRegressionWS1(c(end_vel, end_pos), i)
        
        c_next = end_vel - 3*cur_coef[5]*(time_splits[4]-time_splits[3])^2 - 2*cur_coef[6]*(time_splits[4]-time_splits[3])
        d_next = end_pos - cur_coef[5]*(time_splits[4]-time_splits[3])^3 - cur_coef[6]*(time_splits[4]-time_splits[3])^2 - c_next*(time_splits[4]-time_splits[3])
        cur_path_next = rbind(cur_path_next, c(unname(cur_coef[5:6]), c_next, d_next))
        
        c_second = cur_path_next[i,3] - 3*cur_coef[3]*(time_splits[3]-time_splits[2])^2 - 2*cur_coef[4]*(time_splits[3]-time_splits[2])
        d_second = cur_path_next[i,4] - cur_coef[3]*(time_splits[3]-time_splits[2])^3 - cur_coef[4]*(time_splits[3]-time_splits[2])^2 - c_second*(time_splits[3]-time_splits[2])
        cur_path_second = rbind(cur_path_second, c(unname(cur_coef[3:4]), c_second, d_second))
        
        c_first = cur_path_second[i,3] - 3*cur_coef[1]*(time_splits[2]-time_splits[1])^2 - 2*cur_coef[2]*(time_splits[2]-time_splits[1])
        d_first = cur_path_second[i,4] - cur_coef[1]*(time_splits[2]-time_splits[1])^3 - cur_coef[2]*(time_splits[2]-time_splits[1])^2 - c_first*(time_splits[2]-time_splits[1])
        cur_path_first = rbind(cur_path_first, c(unname(cur_coef[1:2]), c_first, d_first))
        
      }
      
      est_path[[1]] = cur_path_first
      est_path[[2]] = cur_path_next
      est_path = append(est_path, list(cur_path_second), after = 1)
      
      cur_est_pos = rbind(evaluateCubic(Data$t[split_labels == 1],Cubic = est_path[[1]], startTime = time_splits[1]),
                          evaluateCubic(Data$t[split_labels == 2],Cubic = est_path[[2]], startTime = time_splits[2]),
                          evaluateCubic(Data$t[split_labels == 3],Cubic = est_path[[3]], startTime = time_splits[3]))
      
      cur_Error_mat = (cur_Data[,-1] - cur_est_pos)^2
      
      Error_mat[split_labels %in% c(1:3),-1] = cur_Error_mat
      
    } else if(worst_split == 2){
      
      cur_Data = Data[split_labels %in% c(1:4),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_path_prev = c()
      cur_path_first = c()
      cur_path_second = c()
      cur_path_next = c()
      
      runRegressionWS2 = function(end_cond, dim){
        
        end_vel = end_cond[1]
        end_pos = end_cond[2]
        
        cur_X = c()
        
        for(j in 1:4){
          cur_a = (cur_Data$t - time_splits[j])^3*(cur_split_labels == j) +
            (2*(time_splits[j+1]-time_splits[j])^3-3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j]))*(cur_split_labels %in% c(1:j))
          
          cur_b = (cur_Data$t - time_splits[j])^2*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^2-2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j]))*(cur_split_labels %in% c(1:j))
          
          cur_X = cbind(cur_X, cur_a, cur_b)
          
        }
        
        cur_Y = cur_Data[,dim+1] - (end_pos + end_vel*(cur_Data$t - time_splits[5]))
        
        cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
        
        cur_coef
        
      }
      
      for(i in 1:n_dim){
        
        end_pos = est_path[[4]][i,4]
        end_vel = est_path[[4]][i,3]
        
        cur_coef = runRegressionWS2(c(end_vel, end_pos), i)
        
        c_next = end_vel - 3*cur_coef[7]*(time_splits[5]-time_splits[4])^2 - 2*cur_coef[8]*(time_splits[5]-time_splits[4])
        d_next = end_pos - cur_coef[7]*(time_splits[5]-time_splits[4])^3 - cur_coef[8]*(time_splits[5]-time_splits[4])^2 - c_next*(time_splits[5]-time_splits[4])
        cur_path_next = rbind(cur_path_next, c(unname(cur_coef[7:8]), c_next, d_next))
        
        c_second = cur_path_next[i,3] - 3*cur_coef[5]*(time_splits[4]-time_splits[3])^2 - 2*cur_coef[6]*(time_splits[4]-time_splits[3])
        d_second = cur_path_next[i,4] - cur_coef[5]*(time_splits[4]-time_splits[3])^3 - cur_coef[6]*(time_splits[4]-time_splits[3])^2 - c_second*(time_splits[4]-time_splits[3])
        cur_path_second = rbind(cur_path_second, c(unname(cur_coef[5:6]), c_second, d_second))
        
        c_first = cur_path_second[i,3] - 3*cur_coef[3]*(time_splits[3]-time_splits[2])^2 - 2*cur_coef[4]*(time_splits[3]-time_splits[2])
        d_first = cur_path_second[i,4] - cur_coef[3]*(time_splits[3]-time_splits[2])^3 - cur_coef[4]*(time_splits[3]-time_splits[2])^2 - c_first*(time_splits[3]-time_splits[2])
        cur_path_first = rbind(cur_path_first, c(unname(cur_coef[3:4]), c_first, d_first))
        
        c_prev = cur_path_first[i,3] - 3*cur_coef[1]*(time_splits[2]-time_splits[1])^2 - 2*cur_coef[2]*(time_splits[2]-time_splits[1])
        d_prev = cur_path_first[i,4] - cur_coef[1]*(time_splits[2]-time_splits[1])^3 - cur_coef[2]*(time_splits[2]-time_splits[1])^2 - c_prev*(time_splits[2]-time_splits[1])
        cur_path_prev = rbind(cur_path_prev, c(unname(cur_coef[1:2]), c_prev, d_prev))
        
      }
      
      est_path[[1]] = cur_path_prev
      est_path[[2]] = cur_path_first
      est_path[[3]] = cur_path_next
      est_path = append(est_path, list(cur_path_second), after = 2)
      
      cur_est_pos = rbind(evaluateCubic(Data$t[split_labels == 1],Cubic = est_path[[1]], startTime = time_splits[1]),
                          evaluateCubic(Data$t[split_labels == 2],Cubic = est_path[[2]], startTime = time_splits[2]),
                          evaluateCubic(Data$t[split_labels == 3],Cubic = est_path[[3]], startTime = time_splits[3]),
                          evaluateCubic(Data$t[split_labels == 4],Cubic = est_path[[4]], startTime = time_splits[4]))
      
      cur_Error_mat = (cur_Data[,-1] - cur_est_pos)^2
      
      Error_mat[split_labels %in% c(1:4),-1] = cur_Error_mat
      
    } else if(worst_split == length(time_splits) - 3){
      
      cur_Data = Data[split_labels %in% c(-1:2 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_path_prev = c()
      cur_path_first = c()
      cur_path_second = c()
      cur_path_next = c()
      
      runRegressionWSL3 = function(init_cond, dim){
        
        init_vel = init_cond[1]
        init_pos = init_cond[2]
        
        cur_X = c()
        
        for(j in -1:2 + worst_split){
          cur_a = (cur_Data$t - time_splits[j])^3*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^3+3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(length(time_splits)-1)))*(j != (length(time_splits)-1))
          
          cur_b = (cur_Data$t - time_splits[j])^2*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^2+2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(length(time_splits)-1)))*(j != (length(time_splits)-1))
          
          cur_X = cbind(cur_X, cur_a, cur_b)
          
        }
        
        cur_Y = cur_Data[,dim+1] - (init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1]))
        
        cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
        
        cur_coef
        
      }
      
      for(i in 1:n_dim){
        
        init_pos = est_path[[worst_split - 1]][i,4]
        init_vel = est_path[[worst_split - 1]][i,3]
        
        cur_coef = runRegressionWSL3(c(init_vel, init_pos), i)
        
        cur_path_prev = rbind(cur_path_prev, c(unname(cur_coef[1:2]), init_vel, init_pos))
        cur_path_first = rbind(cur_path_first, c(unname(cur_coef[3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1])))
        cur_path_second = rbind(cur_path_second, c(unname(cur_coef[5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split])))
        cur_path_next = rbind(cur_path_next, c(unname(cur_coef[7:8]), evaluateCubicVelocity(time_splits[worst_split + 2], cur_path_second[i,], time_splits[worst_split + 1]), evaluateCubic(time_splits[worst_split + 2], cur_path_second[i,], time_splits[worst_split + 1])))
        
      }
      
      est_path[[worst_split - 1]] = cur_path_prev
      est_path[[worst_split]] = cur_path_first
      est_path[[worst_split + 1]] = cur_path_next
      est_path = append(est_path, list(cur_path_second), after = worst_split)
      
      
      cur_est_pos = rbind(evaluateCubic(Data$t[split_labels == worst_split - 1],Cubic = est_path[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                          evaluateCubic(Data$t[split_labels == worst_split],Cubic = est_path[[worst_split]], startTime = time_splits[worst_split]),
                          evaluateCubic(Data$t[split_labels == worst_split + 1],Cubic = est_path[[worst_split + 1]], startTime = time_splits[worst_split + 1]),
                          evaluateCubic(Data$t[split_labels == worst_split + 2],Cubic = est_path[[worst_split + 2]], startTime = time_splits[worst_split + 2]))
      
      cur_Error_mat = (cur_Data[,-1] - cur_est_pos)^2
      
      Error_mat[split_labels %in% c(-1:2 + worst_split),-1] = cur_Error_mat
      
    } else if(worst_split == length(time_splits) - 2){
      
      cur_Data = Data[split_labels %in% c(-1:1 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_path_prev = c()
      cur_path_first = c()
      cur_path_second = c()
      
      runRegressionWSL2 = function(init_cond, dim){
        
        init_vel = init_cond[1]
        init_pos = init_cond[2]
        
        cur_X = c()
        
        for(j in -1:1 + worst_split){
          cur_a = (cur_Data$t - time_splits[j])^3*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^3+3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(length(time_splits)-1)))*(j != (length(time_splits)-1))
          
          cur_b = (cur_Data$t - time_splits[j])^2*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^2+2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(length(time_splits)-1)))*(j != (length(time_splits)-1))
          
          cur_X = cbind(cur_X, cur_a, cur_b)
          
        }
        
        cur_Y = cur_Data[,dim+1] - (init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1]))
        
        cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
        
        cur_coef
        
      }
      
      
      for(i in 1:n_dim){
        
        init_pos = est_path[[worst_split - 1]][i,4]
        init_vel = est_path[[worst_split - 1]][i,3]
        
        cur_coef = runRegressionWSL2(c(init_vel, init_pos), i)
        
        cur_path_prev = rbind(cur_path_prev, c(unname(cur_coef[1:2]), init_vel, init_pos))
        cur_path_first = rbind(cur_path_first, c(unname(cur_coef[3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1])))
        cur_path_second = rbind(cur_path_second, c(unname(cur_coef[5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split])))
        
      }
      
      
      est_path[[worst_split - 1]] = cur_path_prev
      est_path[[worst_split]] = cur_path_first
      est_path = append(est_path, list(cur_path_second), after = worst_split)
      
      cur_est_pos = rbind(evaluateCubic(Data$t[split_labels == worst_split - 1],Cubic = est_path[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                          evaluateCubic(Data$t[split_labels == worst_split],Cubic = est_path[[worst_split]], startTime = time_splits[worst_split]),
                          evaluateCubic(Data$t[split_labels == worst_split + 1],Cubic = est_path[[worst_split + 1]], startTime = time_splits[worst_split + 1]))
      
      cur_Error_mat = (cur_Data[,-1] - cur_est_pos)^2
      
      Error_mat[split_labels %in% c(-1:1 + worst_split),-1] = cur_Error_mat
      
    } else{
      
      cur_Data = Data[split_labels %in% c(-1:2 + worst_split),]
      cur_split_labels = cut(x = cur_Data$t, breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
      cur_path_prev = c()
      cur_path_first = c()
      cur_path_second = c()
      cur_path_next = c()
      
      runRegressionGen = function(cond, dim){
        
        init_vel = cond[1]
        init_pos = cond[2]
        end_vel = cond[3]
        end_pos = cond[4]
        
        cur_X = c()
        
        for(j in -1:1 + worst_split){
          cur_a = (cur_Data$t - time_splits[j])^3*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^3+3*(time_splits[j+1]-time_splits[j])^2*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(worst_split + 2))) +
            (((2*(time_splits[j+1] - time_splits[j])^3 + 6*(time_splits[j+1] - time_splits[j])^2*(time_splits[worst_split+2] - time_splits[j+1]) + 3*(time_splits[j+1] - time_splits[j])^2*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^3))*(cur_Data$t - time_splits[worst_split+2])^3 -
               ((3*(time_splits[j+1] - time_splits[j])^3 + 9*(time_splits[j+1] - time_splits[j])^2*(time_splits[worst_split+2] - time_splits[j+1]) + 6*(time_splits[j+1] - time_splits[j])^2*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^2))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2)
          
          cur_b = (cur_Data$t - time_splits[j])^2*(cur_split_labels == j) +
            ((time_splits[j+1]-time_splits[j])^2+2*(time_splits[j+1]-time_splits[j])*(cur_Data$t - time_splits[j+1]))*(cur_split_labels %in% c((j+1):(worst_split + 2))) +
            (((2*(time_splits[j+1] - time_splits[j])^2 + 4*(time_splits[j+1] - time_splits[j])*(time_splits[worst_split+2] - time_splits[j+1]) + 2*(time_splits[j+1] - time_splits[j])*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^3))*(cur_Data$t - time_splits[worst_split+2])^3 -
               ((3*(time_splits[j+1] - time_splits[j])^2 + 6*(time_splits[j+1] - time_splits[j])*(time_splits[worst_split+2] - time_splits[j+1]) + 4*(time_splits[j+1] - time_splits[j])*(time_splits[worst_split+3] - time_splits[worst_split+2]))/((time_splits[worst_split+3] - time_splits[worst_split+2])^2))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2)
          
          
          cur_X = cbind(cur_X, cur_a, cur_b)
          
        }
        
        cur_Y = cur_Data[,dim+1] - (init_pos + init_vel*(cur_Data$t - time_splits[worst_split-1]) +
                                      (((end_vel + init_vel)/(time_splits[worst_split+3] - time_splits[worst_split+2])^2 + (2*init_vel*(time_splits[worst_split+2] - time_splits[worst_split-1]) - 2*(end_pos-init_pos))/(time_splits[worst_split+3] - time_splits[worst_split+2])^3)*(cur_Data$t - time_splits[worst_split+2])^3 +
                                         ((3*(end_pos - init_pos) - 3*init_vel*(time_splits[worst_split+2] - time_splits[worst_split-1]))/(time_splits[worst_split+3] - time_splits[worst_split+2])^2 - (end_vel + 2*init_vel)/(time_splits[worst_split+3] - time_splits[worst_split+2]))*(cur_Data$t - time_splits[worst_split+2])^2)*(cur_split_labels == worst_split + 2))
        
        cur_coef = unname(RcppArmadillo::fastLmPure(cur_X, cur_Y)$coefficients)
        
        cur_coef
        
      }
      
      for(i in 1:n_dim){
        
        init_pos = est_path[[worst_split - 1]][i,4]
        init_vel = est_path[[worst_split - 1]][i,3]
        end_pos = est_path[[worst_split + 2]][i,4]
        end_vel = est_path[[worst_split + 2]][i,3]
        
        cur_coef = runRegressionGen(c(init_vel, init_pos, end_vel, end_pos), i)
        
        cur_path_prev = rbind(cur_path_prev, c(unname(cur_coef[1:2]), init_vel, init_pos))
        cur_path_first = rbind(cur_path_first, c(unname(cur_coef[3:4]), evaluateCubicVelocity(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1]), evaluateCubic(time_splits[worst_split], cur_path_prev[i,], time_splits[worst_split - 1])))
        cur_path_second = rbind(cur_path_second, c(unname(cur_coef[5:6]), evaluateCubicVelocity(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split]), evaluateCubic(time_splits[worst_split + 1], cur_path_first[i,], time_splits[worst_split])))
        
        init_vel_next = evaluateCubicVelocity(time_splits[worst_split + 2], cur_path_second[i,], time_splits[worst_split+1])
        init_pos_next = evaluateCubic(time_splits[worst_split + 2], cur_path_second[i,], time_splits[worst_split+1])
        diffTime_next = time_splits[worst_split+3] - time_splits[worst_split+2]
        
        cur_path_next = rbind(cur_path_next, c((init_vel_next + end_vel)/(diffTime_next^2) - 2*(end_pos - init_pos_next)/(diffTime_next^3),
                                               (end_vel - init_vel_next)/(2*diffTime_next) - 3*(init_vel_next + end_vel)/(2*diffTime_next) + 3*(end_pos - init_pos_next)/(diffTime_next^2),
                                               init_vel_next, init_pos_next))
        
      }
      
      est_path[[worst_split - 1]] = cur_path_prev
      est_path[[worst_split]] = cur_path_first
      est_path[[worst_split + 1]] = cur_path_next
      est_path = append(est_path, list(cur_path_second), after = worst_split)
      
      cur_est_pos = rbind(evaluateCubic(Data$t[split_labels == worst_split - 1],Cubic = est_path[[worst_split - 1]], startTime = time_splits[worst_split - 1]),
                          evaluateCubic(Data$t[split_labels == worst_split],Cubic = est_path[[worst_split]], startTime = time_splits[worst_split]),
                          evaluateCubic(Data$t[split_labels == worst_split + 1],Cubic = est_path[[worst_split + 1]], startTime = time_splits[worst_split + 1]),
                          evaluateCubic(Data$t[split_labels == worst_split + 2],Cubic = est_path[[worst_split + 2]], startTime = time_splits[worst_split + 2]))
      
      cur_Error_mat = (cur_Data[,-1] - cur_est_pos)^2
      
      Error_mat[split_labels %in% c(-1:2 + worst_split),-1] = cur_Error_mat
      
    }
    
    pos_splits = unname(table(split_labels_unique) >= 8)
    
    all_splits_est_path[[length(all_splits_est_path)+1]] = list(Path = est_path, TimeSplits = time_splits, Error = Error_mat, n_splits = length(est_path) - 1, PosSplits = pos_splits)
    
    svMisc::progress(length(est_path) - 1, max.value = nsplits, progress.bar = T, console = T)
    
    if(sum(pos_splits) == 0){
      break
    }
    
  }
  
  getBestSplitsPath(all_splits_est_path, show = F, V_Smooth = F, type = "BIC")
  
}
