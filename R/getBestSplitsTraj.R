#' getLikelihoodTraj
#'
#' @param EstTraj an EstimatedTraj object.
#'
#' @returns a scalar.
#' @export
#' 
getLikelihoodTraj = function(EstTraj){
  
  n_models = nrow(EstTraj$Traj[[1]])
  
  SSE_traj_Func = function(i){
    EndPosSSE = c((evaluateCubic(EstTraj$TimeSplits[i+1], EstTraj$Traj[[i]], EstTraj$TimeSplits[i]) - EstTraj$Traj[[i]][,4])^2)
    EndVelSSE = c(evaluateCubicVelocity(EstTraj$TimeSplits[i+1], EstTraj$Traj[[i]], EstTraj$TimeSplits[i])^2)
    
    return(c(EndPosSSE, EndVelSSE))
  }
  
  SSE_traj = sapply(1:length(EstTraj$Traj), FUN = SSE_traj_Func)
  
  est_k2_traj = mean(SSE_traj)
  
  getLike = function(i){
    
    likelihood = c()
    timeDiff = EstTraj$TimeSplits[i+1] - EstTraj$TimeSplits[i]
    
    if(i == 1){
      
      cur_Sigma = (2*est_k2_traj)*matrix(c((timeDiff^2+2)/(timeDiff^6)        , (-3*timeDiff^2 - 6)/(2*timeDiff^5), (1)/(2*timeDiff^2),
                                           (-3*timeDiff^2 - 6)/(2*timeDiff^5) , (5*timeDiff^2 + 9)/(2*timeDiff^4) , (-1)/(timeDiff),
                                           (1)/(2*timeDiff^2)                    , (-1)/(timeDiff)                      , 1/2), 
                                         nrow = 3, byrow = T)
      
      for(j in 1:n_models){
        
        likelihood = c(likelihood, mvtnorm::dmvnorm(EstTraj$Traj[[i]][j,1:3], mean = c(0,0,0), sigma = cur_Sigma, log = T, checkSymmetry = F))
        
      }
      
    } else{
      
      cur_Sigma = (2*est_k2_traj)*matrix(c((timeDiff^2+4)/(2*timeDiff^6) , (-timeDiff^2-6)/(2*timeDiff^5),
                                           (-timeDiff^2-6)/(2*timeDiff^5), (timeDiff^2+9)/(2*timeDiff^4)), 
                                         nrow = 2, byrow = T)
      
      for(j in 1:n_models){
        
        likelihood = c(likelihood, mvtnorm::dmvnorm(EstTraj$Traj[[i]][j,1:2], mean = EstTraj$Traj[[i]][j,3]*c(1/timeDiff^2, -2/timeDiff), sigma = cur_Sigma, log = T, checkSymmetry = F))
        
      }
      
    }
    
    return(sum(likelihood))
    
  }
  
  sum(sapply(1:length(EstTraj$Traj), FUN = getLike))
  
}


#' getBestSplitsTraj
#'
#' @param EstTrajList a list of EstimatedTraj objects. 
#' They should all be related using the same path with an increasing number of splits.
#' @param EstPath an EstimatedPath object. This is the path all the trajectories where estimated from.
#' @param show a boolean. TRUE if the criterion values should be displayed.
#' @param type a string. The model criterion. 
#' Must pick from the following: BIC, BICOpt, BICPos, AIC, VarPos. 
#' BIC and AIC can be used with either trajectory method. BICOpt, BICPos, and VarPos can only be used with Optimization.
#' @param V_Smooth a boolean. TRUE if the trajectory cubic splines had C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedTraj object.
#' @export
#'
getBestSplitsTraj = function(EstTrajList, EstPath, show = F, type = "BIC", V_Smooth = T){
  
  if(type == "BIC"){
    
    getBIC = function(cur_traj){
      
      n_models = nrow(cur_traj$Traj[[1]])
      n_points = nrow(cur_traj$Error)*(ncol(cur_traj$Error)-1)
      n_param = n_models*(4 + V_Smooth*(length(cur_traj$Traj)-1) + (!V_Smooth)*2*(length(cur_traj$Traj)-1))
      
      return(n_param*log(n_points) + n_points*log(stats::median(rowSums(cur_traj$Error[,-1]))) - 2*getLikelihoodTraj(EstTraj = cur_traj))
    }
    
    BIC = sapply(EstTrajList, getBIC)
    
    best_index = which.min(BIC)
    
    print(stringr::str_c("Best Number of Splits was ", best_index-1,"."))
    if(show){
      print(BIC)
    }
    
  }
  
  if(type == "BICOpt"){
    
    getBIC = function(cur_traj){
      
      n_models = nrow(cur_traj$Traj[[1]])
      n_points = nrow(EstPath$Error)*(ncol(EstPath$Error)-1)
      n_param = n_models*(4 + V_Smooth*(length(cur_traj$Traj)-1) + (!V_Smooth)*2*(length(cur_traj$Traj)-1))
      
      return(n_param*log(n_points) + n_points*log(sum(cur_traj$Error)/(max(cur_traj$TimeSplits)-min(cur_traj$TimeSplits))) - 2*getLikelihoodTraj(EstTraj = cur_traj))

    }
    
    
    BIC = sapply(EstTrajList, getBIC)
    
    best_index = which.min(BIC)
    
    print(stringr::str_c("Best Number of Splits was ", best_index+2,"."))
    if(show){
      print(BIC)
    }
    
  }
  
  if(type == "BICPos"){
    
    getBICPos = function(cur_traj){
      
      t_seq = cur_traj$Error[,1]
      n_points = length(t_seq)
      n_models = nrow(cur_traj$Traj[[1]])
      n_param = n_models*(4 + V_Smooth*(length(cur_traj)-1) + (!V_Smooth)*2*(length(cur_traj)-1))
      TimeSplits = cur_traj$TimeSplits
      Cov = cur_traj$Cov
      Sigma2 = cur_traj$Sigma2
      EstVelCov = cur_traj$EstVelCovMat
      n_dim = nrow(EstVelCov)/n_points
      Error = cur_traj$Error
      
      evaluateIndividualTraj = function(t){
        cur_traj_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_traj_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
        
        (mean(apply(matrix(1:n_models), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_traj_index]][1:4 + (x-1)*4, 1:4 + (x-1)*4] %*% trans_mat)[1,1] + sum(EstVelCov[(0:(n_dim-1))*n_points + which(t_seq == t), (0:(n_dim-1))*n_points + which(t_seq == t)]) + Sigma2[[cur_traj_index]][x]})))*(mean(colMeans(Error[,-1])))
        
      }
      
      cur_PosVar = unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualTraj), probs = c(0.95)))
      
      cur_BIC = (n_param)*log(n_points) + n_points*log(stats::median(unlist(as.vector(cur_traj$Error[,-1]))))
      
      cur_BIC / cur_PosVar
      
      
    }
    
    BIC_Pos = sapply(EstTrajList, getBICPos)
    
    best_index = which.min(BIC_Pos)
    
    print(stringr::str_c("Best Number of Splits was ", best_index,"."))
    if(show){
      print(BIC_Pos)
    }
    
  }
  
  if(type == "AIC"){
    
    getAIC = function(cur_traj){
      
      n_models = nrow(cur_traj$Traj[[1]])
      n_points = nrow(cur_traj$Error)*(ncol(cur_traj$Error)-1)
      
      return(2*(2*n_models*(length(cur_traj$Traj)+1)) + n_points*log(stats::median(unlist(as.vector(cur_traj$Error[,1])))) - 2*getLikelihoodTraj(EstTraj = cur_traj))
      
    }
    
    AIC = sapply(EstTrajList, getAIC)
    
    best_index = which.min(AIC)
    
    print(stringr::str_c("Best Number of Splits was ", best_index-1,"."))
    if(show){
      print(AIC)
    }
    
  }
  
  if(type == "VarPos"){
    
    getMaxEstPosVar = function(cur_traj){
      
      t_seq = seq(min(cur_traj$TimeSplits), max(cur_traj$TimeSplits), by = 0.01)
      n_points = length(t_seq)
      n_models = nrow(cur_traj$Traj[[1]])
      TimeSplits = cur_traj$TimeSplits
      Cov = cur_traj$Cov
      Sigma2 = cur_traj$Sigma2
      
      evaluateIndividualTraj = function(t){
        cur_traj_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_traj_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
        
        (mean(apply(matrix(1:n_models), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_traj_index]][1:4 + (x-1)*4, 1:4 + (x-1)*4] %*% trans_mat)[1,1] + Sigma2[[cur_traj_index]]})))
        
      }
      
      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualTraj), probs = c(0.95)))
      
    }
    
    max_PosVar = sapply(EstTrajList, FUN = getMaxEstPosVar, simplify = T)
    
    best_index = which.min(max_PosVar)
    
    print(stringr::str_c("Best Number of Splits was ", best_index,"."))
    if(show){
      print(max_PosVar)
    }
    
  }
  
  EstTrajList[[best_index]]
  
}

