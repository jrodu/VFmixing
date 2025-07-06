#' 4th-Order Runge-Kutta Method
#'
#' @param startTime a numeric scalar.
#' @param startPos a numeric vector.
#' @param baseVectorFields a baseVectorField function.
#' @param TrajList a list of numeric matrices. The list of coefficient matrices of the EstimatedTrajectory object.
#' @param TrajTimeSplits a numeric vector. The nodes of the EstimatedTrajectory object in addition with the starting and ending time.
#' @param endTime a numeric scalar.
#' @param t_step a positive numeric scalar. The time grid size for the Runge-Kutta method. 
#'
#' @returns a numeric matrix. Each row is a time and position of the path through the trajectory-weighted vector field.
#' @export
#' 
#' @examples
#' RungeKutta(startTime = 0, startPos = c(1,1), 
#'            baseVectorFields = baseVectorFields,
#'            TrajList = ExTraj$Traj, TrajTimeSplits = ExTraj$TimeSplits, 
#'            endTime = 5, t_step = 0.01)
#'
RungeKutta = function(startTime, startPos, baseVectorFields, TrajList, TrajTimeSplits, endTime, t_step = 0.01){
  
  n_dim = length(startPos)
  
  t_sim = seq(startTime, endTime, t_step)
  
  pos_sim = matrix(startPos, nrow = 1, byrow = T)
  i=1
  for(i in 1:(length(t_sim)-1)){
    
    cur_t = t_sim[i]
    cur_pos = pos_sim[i,]
    
    k1 = TrajWeightedBaseVectorFields(cur_t, cur_pos, baseVectorFields, TrajList, TrajTimeSplits)
    k2 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k1, baseVectorFields, TrajList, TrajTimeSplits)
    k3 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k2, baseVectorFields, TrajList, TrajTimeSplits)
    k4 = TrajWeightedBaseVectorFields(cur_t + t_step, cur_pos + t_step*k3, baseVectorFields, TrajList, TrajTimeSplits)
    
    pos_sim = rbind(pos_sim, cur_pos + t_step*(k1+2*k2+2*k3+k4)/6)
    
  }
  
  full_sim = data.frame(cbind(t_sim, pos_sim))
  
  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))
  
  full_sim
}