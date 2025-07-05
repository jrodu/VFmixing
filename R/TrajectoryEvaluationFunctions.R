#' getEstTrajValue
#'
#' @param t a numeric vector. The times at which the trajectory's value should be evaluated.
#' @param EstTraj an EstimatedTrajectory object.
#'
#' @returns a numeric matrix. Each row corresponds to an evaluated value.
#' @export
#'
getEstTrajValue = function(t, EstTraj){
  
  n_models = nrow(EstTraj$Traj[[1]])
  TimeSplits = EstTraj$TimeSplits
  
  IndTime = function(t){
    
    cur_traj_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_traj_index]
    trans_mat = matrix(c(cur_difft^3, cur_difft^2, cur_difft, 1))
    
    cur_traj = EstTraj$Traj[[cur_traj_index]] %*% trans_mat
    
    cur_traj
    
    
  }
  
  matrix(t(sapply(X = t, FUN = IndTime)), ncol = n_models)
  
}

#' getEstTrajSlope
#'
#' @param t a numeric vector. The times at which the trajectory's slope should be evaluated.
#' @param EstTraj an EstimatedTrajectory object.
#'
#' @returns a numeric matrix. Each row corresponds to an evaluated slope.
#' @export
#'
getEstTrajSlope = function(t, EstTraj){
  
  n_models = nrow(EstTraj$Traj[[1]])
  TimeSplits = EstTraj$TimeSplits
  
  IndTime = function(t){
    
    cur_traj_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_difft = t - TimeSplits[cur_traj_index]
    trans_mat = matrix(c(3*cur_difft^2, 2*cur_difft, 1, 0))
    
    cur_traj = EstTraj$Traj[[cur_traj_index]] %*% trans_mat
    
    cur_traj
    
  }
  
  matrix(t(sapply(X = t, FUN = IndTime)), ncol = n_models)
  
}