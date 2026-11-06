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
#' Must pick from the following: BIC, BICOS, AIC, VarPos.
#' BIC and AIC can be used with any trajectory method. "BICOS" is kept for compatibility and is the same as "BIC".
#' The old name "BICOpt" is accepted as an alias for "BICOS".
#' VarPos needs a covariance structure, so it can only be used with OneShot.
#' "BICPos" is no longer available: it needed a velocity covariance matrix that no trajectory method produces.
#' @param V_Smooth a boolean. TRUE if the trajectory cubic splines had C2 smoothness. FALSE if C1.
#'
#' @returns an EstimatedTraj object.
#' @export
#'
getBestSplitsTraj = function(EstTrajList, EstPath, show = F, type = "BIC", V_Smooth = T){

  if(identical(type, "BICOpt")){
    type = "BICOS"
  }
  
  if(identical(type, "BICPos")){
    stop("type = \"BICPos\" is no longer available: it needed a velocity covariance matrix ",
         "that no trajectory method produces. Use \"BIC\" or \"VarPos\" instead.", call. = FALSE)
  }

  valid_types = c("BIC", "BICOS", "AIC", "VarPos")

  if(length(type) != 1 || !(type %in% valid_types)){
    stop("type must be one of: ", paste(valid_types, collapse = ", "), call. = FALSE)
  }
  
  if(type == "VarPos" && is.null(EstTrajList[[1]]$Cov)){
    stop("type = \"VarPos\" needs a trajectory with a covariance structure (OneShot). ",
         "Use \"BIC\" or \"AIC\" instead.", call. = FALSE)
  }

  n_dim = nrow(EstPath$Path[[1]])
  
  # Effective sample size: the number of free parameters in the path the 
  # trajectories were estimated from.
  EDF = n_dim*splineFreeParams(length(EstPath$Path), inferPathVSmooth(EstPath, default = V_Smooth))

  if(type %in% c("BIC", "BICOS")){
    
    getBIC = function(cur_traj){
      
      n_models = nrow(cur_traj$Traj[[1]])
      n_param = n_models*splineFreeParams(length(cur_traj$Traj), V_Smooth)
      
      n_param*log(EDF) + EDF*log(trajErrorStat(cur_traj))
    }
    
    criterion = sapply(EstTrajList, getBIC)
    
  }
  
  if(type == "AIC"){
    
    getAIC = function(cur_traj){
      
      n_models = nrow(cur_traj$Traj[[1]])
      
      2*(2*n_models*(length(cur_traj$Traj)+1)) + EDF*log(trajErrorStat(cur_traj))
      
    }
    
    criterion = sapply(EstTrajList, getAIC)
    
  }
  
  if(type == "VarPos"){
    
    getMaxEstPosVar = function(cur_traj){
      
      t_seq = seq(min(cur_traj$TimeSplits), max(cur_traj$TimeSplits), by = 0.01)
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
    
    criterion = sapply(EstTrajList, FUN = getMaxEstPosVar, simplify = T)
    
  }
  
  best_index = which.min(criterion)
  
  reportBestRegions(length(EstTrajList[[best_index]]$Traj))
  if(show){
    print(criterion)
  }
  
  EstTrajList[[best_index]]
  
}
