#' Path Model Selection
#'
#' @param EstPath a list of EstimatedPath Objects. These should all be related using the same data with an increasing number of splits.
#' @param show a boolean. TRUE if the criterion values should be displayed.
#' @param V_Smooth a boolean. TRUE if the path cubic spline had C2 smoothness. FALSE if C1.
#' @param type a string. The model criterion. 
#' Must pick from the following: BIC, BICVel, AIC, VarPos, VarVel, VarPosVel. 
#' BIC and AIC can be used with any path method. BICVel, VarPos, VarVel, VarPosVel need a covariance structure,
#' so they can only be used with FullRegression, Projection and FullFours.
#'
#' @returns an EstimatedPath object.
#' @export
#'
getBestSplitsPath = function(EstPath, show = F, V_Smooth = T, type = "BIC"){

  valid_types = c("BIC", "BICVel", "AIC", "VarPos", "VarVel", "VarPosVel")

  if(length(type) != 1 || !(type %in% valid_types)){
    stop("type must be one of: ", paste(valid_types, collapse = ", "), call. = FALSE)
  }

  if(type %in% c("BICVel", "VarPos", "VarVel", "VarPosVel") && is.null(EstPath[[1]]$Cov)){
    stop("type = \"", type, "\" needs a path with a covariance structure ",
         "(FullRegression, Projection or FullFours). Use \"BIC\" or \"AIC\" instead.", call. = FALSE)
  }

  if(type == "BIC"){
    
    getBIC = function(cur_path){
      
      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)
      n_param = n_dim*splineFreeParams(length(cur_path$Path), V_Smooth)
      
      return((n_param)*log(n_points) + n_points*log(stats::median(rowSums(cur_path$Error[,-1]))))
      
    }
    
    BIC = sapply(EstPath, getBIC)
    
    best_index = which.min(BIC)
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(BIC)
    }
    
  }
  
  if(type == "BICVel"){
    
    getBICVel = function(cur_path){
      
      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      Error = cur_path$Error
      n_param = n_dim*splineFreeParams(length(cur_path$Path), V_Smooth)

      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)
        
        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))
        
      }
      
      cur_VelVar = unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))
      
      cur_BIC = (n_param)*log(n_points) + n_points*log(stats::median(unlist(as.vector(cur_path$Error[,-1]))))
      
      c(cur_BIC, cur_VelVar)
      
    }
    
    BIC_Vel = sapply(EstPath, getBICVel)
    
    # Only consider models before the velocity variance first jumps by more
    # than a factor of e. If it never jumps, consider every model.
    vel_cutoff = which(diff(log(BIC_Vel[2,])) > 1)[1]

    if(is.na(vel_cutoff)){
      vel_cutoff = ncol(BIC_Vel)
    }

    best_index = which.min(BIC_Vel[1, 1:vel_cutoff])
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(BIC_Vel)
    }
    
    
  }
  
  if(type == "AIC"){
    
    getAIC = function(cur_path){
      
      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)
      
      return(2*(2*n_dim*(length(cur_path$Path)+1)) + n_points*log(stats::median(unlist(as.vector(cur_path$Error[-1])))))
    }
    
    AIC = sapply(EstPath, getAIC)
    
    best_index = which.min(AIC)
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(AIC)
    }
    
  }
  
  if(type == "VarPos"){
    
    getMaxEstPosVar = function(cur_path){
      
      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      Error = cur_path$Error
      
      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
        
        (mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]})))
        
      }
      
      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))
      
    }
    
    max_PosVar = sapply(EstPath, FUN = getMaxEstPosVar, simplify = T)
    
    best_index = which.min(max_PosVar)
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(max_PosVar)
    }
    
  }
  
  if(type == "VarVel"){
    
    getMaxEstVelVar = function(cur_path){
      
      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      
      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)
        
        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))
        
      }
      
      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))
      
    }
    
    max_VelVar = sapply(EstPath, FUN = getMaxEstVelVar, simplify = T)
    
    best_index = which.min(max_VelVar)
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(max_VelVar)
    }
    
  }
  
  if(type == "VarPosVel"){
    
    getMaxEstPosVelVar = function(cur_path){
      
      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      
      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1) + c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)
        
        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))
        
      }
      
      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))
      
    }
    
    max_PosVelVar = sapply(EstPath, FUN = getMaxEstPosVelVar, simplify = T)
    
    best_index = which.min(max_PosVelVar)
    n_regions = length(EstPath[[best_index]]$Path)
    
    reportBestRegions(n_regions)
    if(show){
      print(max_PosVelVar)
    }
    
  }
  
  EstPath[[best_index]]
  
  
}
