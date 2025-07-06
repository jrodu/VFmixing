#' Generate a random position from an Estimated Path
#' 
#' @description
#' Generating random positions using the normal distrubutions induced by the linear regression.
#' 
#'
#' @param t a numeric vector. The times at which random positions should be generated.
#' @param EstPath an EstimatedPath object. Must have been estimated using the FullRegression or Projection methods.
#' @param random a boolean. If TRUE the positions are random. 
#'
#' @returns a numeric matrix. Each row corresponds to a random position at a time in the t vector.
#' @export
#' 
#' @examples
#' getRandomEstPathPosition(c(0,0.5,1), ExEstPath_FullRegression)
#' 
#'
getRandomEstPathPosition = function(t, EstPath, random = T){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  Cov = EstPath$Cov
  Sigma2 = EstPath$Sigma2
  
  sampleIndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(cur_difft^3, cur_difft^2, cur_difft, 1))
    
    samp_pos = c()
    for(i in 1:n_dim){
      
      est_pos = EstPath$Path[[cur_path_index]][i,] %*% trans_mat
      cur_var = t(trans_mat) %*% Cov[[cur_path_index]][[i]] %*% trans_mat + Sigma2[[cur_path_index]][[i]]
      
      samp_pos = c(samp_pos, stats::rnorm(1, mean = est_pos, sd = sqrt(cur_var)*random))
      
    }
    
    samp_pos
    
  }
  
  matrix(t(sapply(X = t, FUN = sampleIndTime)), ncol = n_dim)
}

#' Generate a random velocity from an Estimated Path
#'
#' @param t a numeric vector. The times at which random velocities should be generated.
#' @param EstPath an EstimatedPath object. Must have been estimated using the FullRegression or Projection methods.
#' @param random a boolean. If TRUE the velocities are random. 
#'
#' @returns a numeric matrix. Each row corresponds to a random velocity at a time in the t vector.
#' @export
#' 
#' @examples
#' getRandomEstPathVelocity(c(0,0.5,1), ExEstPath_FullRegression)
#'
getRandomEstPathVelocity = function(t, EstPath, random = T){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  Cov = EstPath$Cov
  Sigma2 = EstPath$Sigma2
  
  sampleIndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(3*cur_difft^2, 2*cur_difft, 1, 0))
    
    samp_pos = c()
    for(i in 1:n_dim){
      
      est_pos = EstPath$Path[[cur_path_index]][i,] %*% trans_mat
      cur_var = t(trans_mat) %*% Cov[[cur_path_index]][[i]] %*% trans_mat  + Sigma2[[cur_path_index]][[i]]
      
      samp_pos = c(samp_pos, stats::rnorm(1, mean = est_pos, sd = sqrt(cur_var)*random))
      
    }
    
    samp_pos
    
  }
  
  matrix(t(sapply(X = t, FUN = sampleIndTime)), ncol = n_dim)
  
}

#' Generate both a random position and velocity from an Estimated Path
#'
#' @param t a numeric vector. The times at which random positions and velocities should be generated.
#' @param EstPath an EstimatedPath object. Must have been estimated using the FullRegression or Projection methods.
#' @param random a boolean. If TRUE the positions and velocities are random. 
#'
#' @returns a numeric matrix. Each row corresponds to a random position and velocity at a time in the t vector.
#' @export
#' 
#' @examples
#' getRandomEstPathPositionVelocity(c(0,0.5,1), ExEstPath_FullRegression)
#' 
getRandomEstPathPositionVelocity = function(t, EstPath, random = T){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  Cov = EstPath$Cov
  Sigma2 = EstPath$Sigma2
  
  sampleIndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(cur_difft^3, cur_difft^2, cur_difft, 1,
                         3*cur_difft^2, 2*cur_difft, 1, 0), nrow = 2)
    
    samp_pos = c()
    samp_vel = c()
    
    for(i in 1:n_dim){
      
      est_posvel = EstPath$Path[[cur_path_index]][i,] %*% t(trans_mat)
      
      cur_cov = trans_mat %*% Cov[[cur_path_index]][[i]] %*% t(trans_mat) + Sigma2[[cur_path_index]][[i]]*diag(2)
      
      cur_samp = MASS::mvrnorm(n = 1, mu = est_posvel, Sigma = cur_cov*random)
      
      samp_pos = c(samp_pos, cur_samp[1])
      samp_vel = c(samp_vel, cur_samp[2])
      
    }
    
    c(samp_pos, samp_vel)
    
  }
  
  matrix(t(sapply(X = t, FUN = sampleIndTime)), ncol = 2*n_dim)
  
}