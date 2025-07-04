#' Evaluate the position of a path
#'
#' @param t a numeric vector. The times at which the path's position should be evaluated.
#' @param EstPath an EstimatedPath object.
#'
#' @returns a numeric matrix. Each row corresponds to an evaluated position.
#' @export
#'
getEstPathPosition = function(t, EstPath){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  
  IndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(cur_difft^3, cur_difft^2, cur_difft, 1))
    
    cur_pos = EstPath$Path[[cur_path_index]] %*% trans_mat
    
    cur_pos
    
  }
  
  matrix(t(sapply(X = t, FUN = IndTime)), ncol = n_dim)
  
  
}

#' Evaluate the velocity of a path
#'
#' @param t a numeric vector. The times at which the path's velocity should be evaluated.
#' @param EstPath an EstimatedPath object.
#'
#' @returns a numeric matrix. Each row corresponds to an evaluated velocity.
#' @export
#'
getEstPathVelocity = function(t, EstPath){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  
  IndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(3*cur_difft^2, 2*cur_difft, 1, 0))
    
    cur_pos = EstPath$Path[[cur_path_index]] %*% trans_mat
    
    cur_pos
    
  }
  
  matrix(t(sapply(X = t, FUN = IndTime)), ncol = n_dim)
  
  
}

#' Evaluate the acceleration of a path
#'
#' @param t a numeric vector. The times at which the path's acceleration should be evaluated.
#' @param EstPath an EstimatedPath object.
#'
#' @returns a numeric matrix. Each row corresponds to an evaluated acceleration.
#' @export
#'
getEstPathAcceleration = function(t, EstPath){
  
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits = EstPath$TimeSplits
  
  IndTime = function(t){
    
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(6*cur_difft, 2, 0, 0))
    
    cur_pos = EstPath$Path[[cur_path_index]] %*% trans_mat
    
    cur_pos
    
  }
  
  matrix(t(sapply(X = t, FUN = IndTime)), ncol = n_dim)
  
  
}