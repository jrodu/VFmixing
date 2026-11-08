#' Get Estimated Path using the By Fours Method
#' 
#' @description
#' One of the functions to estimate a path. ByFours method only fits the previous segment, the two splitting segments, and the next segment during each iteration. 
#' 
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param nsplits a numeric scalar. The maximum number of splits the spline can split into.
#' @param Model_Selection_Type a string. The model criterion, "BIC" or "AIC". 
#' The other criteria need a covariance structure, which ByFours does not produce.
#' @param return_all a boolean. If TRUE, return the list of every candidate path from the node search instead of the selected one. 
#' Used by [getEstimatedPaths_FullFours()].
#' @param patience NULL or a positive integer. If NULL (the default), every candidate up to `nsplits` is fit.
#' Otherwise, stop the node search after this many consecutive candidates that do not improve the selection criterion.
#' This is faster, but it can miss a better model that only appears after a long run without improvement.
#' @param undersmooth a boolean. If TRUE (the default), return the candidate with twice as many regions as the one
#' `Model_Selection_Type` selects (the first candidate with at least that many, or the last candidate if the search
#' ends first), instead of the selected one. This undersmooths the fit, which reduces the bias of the estimated
#' velocity so that the intervals of the trajectories estimated from it are centred closer to the truth.
#' With `patience`, the criterion's choice is fixed when patience runs out, and the search then continues until
#' it reaches twice that many regions (or returns at once if it already has).
#'
#' @returns an EstimatedPath object.
#' @export
#'
#' @examples
#' getEstimatedPaths_ByFours(ExData, 100)
#'
getEstimatedPaths_ByFours = function(Data, nsplits, Model_Selection_Type = "BIC", return_all = F, patience = NULL, undersmooth = T){

  if(!(Model_Selection_Type %in% c("BIC", "AIC"))){
    stop("Model_Selection_Type must be \"BIC\" or \"AIC\" for the ByFours method.", call. = FALSE)
  }

  checkPatience(patience)

  # ByFours candidates are always C1.
  tracker = pathPatienceTracker(patience, Model_Selection_Type, V_Smooth = F, undersmooth = undersmooth)

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))
  
  # Start from a single region. The local "by fours" update needs at least
  # four regions to work with, so the first candidates are fit with the full
  # C1 spline regression (the same spline space ByFours fits) and split the
  # same way as below: the region with the largest squared error, at the
  # median time of its data.
  time_splits = range(Data$t)
  
  all_splits_est_path = list()
  
  asByFours = function(fit){
    list(Path = fit$Path, TimeSplits = fit$TimeSplits, Error = fit$Error, n_splits = length(fit$Path), PosSplits = fit$PosSplits)
  }
  
  fit = runEstPathStep_FullRegression(Data, time_splits, V_Smooth = F)
  all_splits_est_path[[1]] = asByFours(fit)
  stop_early = !is.null(patience) && tracker$add(all_splits_est_path[[1]])
  
  svMisc::progress(0, max.value = nsplits, progress.bar = T, console = T)
  
  while(!stop_early && length(fit$Path) < 4 && length(time_splits) - 1 <= nsplits){
    
    last = all_splits_est_path[[length(all_splits_est_path)]]
    
    if(sum(last$PosSplits) == 0){
      break
    }
    
    worst_split = which.max(unname(by(rowSums(last$Error[,-1]), fit$SplitLabels, sum))*last$PosSplits)
    new_splits = c(time_splits[1:worst_split], stats::median(Data$t[fit$SplitLabels == worst_split]), time_splits[(worst_split+1):length(time_splits)])
    
    new_fit = tryCatch(runEstPathStep_FullRegression(Data, new_splits, V_Smooth = F), error = function(e){NULL})
    
    if(is.null(new_fit)){
      # This region can't be split; try another one
      all_splits_est_path[[length(all_splits_est_path)]]$PosSplits[worst_split] = FALSE
      next
    }
    
    fit = new_fit
    time_splits = new_splits
    all_splits_est_path[[length(all_splits_est_path) + 1]] = asByFours(fit)
    stop_early = !is.null(patience) && tracker$add(all_splits_est_path[[length(all_splits_est_path)]])
    
    svMisc::progress(length(time_splits) - 2, max.value = nsplits, progress.bar = T, console = T)
    
  }
  
  # Continue with the by-fours updates from the four-region fit
  est_path = fit$Path
  Error_mat = fit$Error
  split_labels = fit$SplitLabels
  split_labels_unique = cut(x = unique(Data$t), breaks = time_splits, include.lowest = T, right = T, labels = 1:(length(time_splits)-1))
  pos_splits = all_splits_est_path[[length(all_splits_est_path)]]$PosSplits
  
  if(length(fit$Path) < 4 || sum(pos_splits) == 0){
    stop_early = TRUE
  }
  

  
  while(!stop_early && length(time_splits)-1 <= nsplits){
    
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
    
    all_splits_est_path[[length(all_splits_est_path)+1]] = list(Path = est_path, TimeSplits = time_splits, Error = Error_mat, n_splits = length(est_path), PosSplits = pos_splits)
    
    svMisc::progress(length(est_path) - 1, max.value = nsplits, progress.bar = T, console = T)

    if(!is.null(patience) && tracker$add(all_splits_est_path[[length(all_splits_est_path)]])){
      break
    }

    if(sum(pos_splits) == 0){
      break
    }

  }

  if(return_all){
    # Keep the criterion's choice at the point patience ran out, for the caller's selection
    attr(all_splits_est_path, "best_index") = if(is.null(patience)) NULL else tracker$best_index()
    return(all_splits_est_path)
  }

  if(!is.null(patience)){
    # The candidates were already scored while searching.
    return(selectBestPath(all_splits_est_path, tracker$values(), Model_Selection_Type,
                          undersmooth = undersmooth, best_index = tracker$best_index()))
  }

  getBestSplitsPath(all_splits_est_path, show = F, V_Smooth = F, type = Model_Selection_Type, undersmooth = undersmooth)

}
