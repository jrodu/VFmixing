#' animatePathTraj
#' 
#' @description
#' A function to animate a particle traveling through a trajectory-weighted vector field.
#' 
#'
#' @param PathList a list of numeric matrices. The list of coefficient matrices of the path.
#' @param TrajList a list of numeric matrices. The list of coefficient matrices of the trajectory. 
#' @param PathTimeSplits a numeric vector. The nodes of the path in addition with the starting and ending time.
#' @param TrajTimeSplits a numeric vector. The nodes of the trajectory in addition with the starting and ending time.
#' @param baseVectorFields a baseVectorFields function.
#' @param n_points the number of points to include in the animation.
#'
#' @returns a ggplot2 object. Can be further turned into an animation.
#' @export
animatePathTraj = function(PathList, TrajList, PathTimeSplits, TrajTimeSplits, baseVectorFields, n_points){
  n_dim = nrow(PathList[[1]])
  n_traj = nrow(TrajList[[1]])
  t_seq = seq(min(PathTimeSplits), max(PathTimeSplits), length.out = n_points)
  
  pos_df = data.frame()
  i=1
  for(i in 1:(length(PathTimeSplits)-1)){
    cur_t_seq = t_seq[t_seq >= PathTimeSplits[i] & t_seq < PathTimeSplits[i+1]]
    
    curPath = PathList[[i]]
    
    cur_pos = evaluateCubic(cur_t_seq, PathList[[i]], startTime = PathTimeSplits[i])
    
    pos_df = rbind(pos_df, cbind(cur_t_seq, cur_pos))
    
  }
  
  names(pos_df) = c('t',stringr::str_c('X', 1:n_dim))
  
  X1_seq = seq(from = floor(min(pos_df$X1)), to = ceiling(max(pos_df$X1)), length.out = 10)
  X2_seq = seq(from = floor(min(pos_df$X2)), to = ceiling(max(pos_df$X2)), length.out = 10)
  velocity_grid = expand.grid(X1_seq, X2_seq)
  names(velocity_grid) = c("X1","X2")
  
  velocity_grid_df = data.frame()
  j=1
  for(j in 1:length(pos_df$t)){
    
    velocity_grid$X1v = 0
    velocity_grid$X2v = 0
    velocity_grid$t = pos_df$t[j]
    i=1
    for(i in 1:nrow(velocity_grid)){
      cur_t = velocity_grid$t[i]
      cur_pos = c(velocity_grid$X1[i], velocity_grid$X2[i])
      
      cur_traj_index = min(c(max(c(sum((TrajTimeSplits - cur_t)<0),1)), length(TrajTimeSplits)-1))
      
      cur_traj = evaluateCubic(t = cur_t, Cubic = TrajList[[cur_traj_index]], startTime = TrajTimeSplits[cur_traj_index])
      
      theorVel = baseVectorFields(t = cur_t, curPos = cur_pos)
      
      theorVel_w = theorVel %*% t(cur_traj)
      
      velocity_grid$X1v[i] = theorVel_w[1]
      velocity_grid$X2v[i] = theorVel_w[2]
    }
    
    velocity_grid_df = rbind(velocity_grid_df, velocity_grid)
    
  }
  
  animation = ggplot2::ggplot() + 
    ggplot2::geom_point(data = pos_df, ggplot2::aes(x = X1, y = X2), size = 1.5) +
    ggplot2::geom_segment(data = velocity_grid_df, ggplot2::aes(x = X1, y = X2, xend = X1 + X1v/5,
                                              yend = X2 + X2v/5), color = 'black',
                 arrow = grid::arrow(length = grid::unit(0.1,"cm")), size = 0.5) + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "aliceblue")) +
    gganimate::transition_time(t) + 
    ggplot2::labs(title = "Current Time: {format(round(frame_time,2),nsmall = 2)}")
  
  animation
  
}