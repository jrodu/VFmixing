#' Visualize Trajectory Functions
#' 
#' @description
#' A function to visualize a series of trajectory functions. 
#' 
#'
#' @param TrajList a list of lists. Each element of the outer list corresponds to a list of cubic matrix coefficients associated with individual cubic functions.
#' @param TimeSplitsList a list of numeric vectors. Each element corresponds to nodes of the cubic spline in the current element of TrajList.
#' @param t_grid_size a non-negative numeric scalar. The time grid size of the plotting.
#'
#' @returns a ggplot2 object. 
#' @export
#'
#' @examples
#' visualizeFullTrajectories(list(ExTraj$Traj),list(ExTraj$TimeSplits), 0.01)
#' 
visualizeFullTrajectories = function(TrajList, TimeSplitsList, t_grid_size){
  Full_Traj_df = data.frame()
  
  for(j in 1:length(TrajList)){
    
    cur_TrajList = TrajList[[j]]
    TimeSplits = TimeSplitsList[[j]]
    
    sub_Traj_df = data.frame()
    
    for(i in 1:length(cur_TrajList)){
      cur_t_seq = seq(TimeSplits[i], TimeSplits[i+1], by = t_grid_size) 
      
      cur_traj_df = data.frame(cbind(rep(cur_t_seq, nrow(cur_TrajList[[i]])), 
                                     c(apply(X = cur_TrajList[[i]], MARGIN = 1, FUN = evaluateCubic, t = cur_t_seq, startTime = TimeSplits[i], simplify = T)),
                                     rep(stringr::str_c("Trajectory", 1:nrow(cur_TrajList[[i]])), each = length(cur_t_seq))))
      names(cur_traj_df) = c("t","Trajectory","Model")
      cur_traj_df$t = as.numeric(cur_traj_df$t)
      cur_traj_df$Trajectory = as.numeric(cur_traj_df$Trajectory)
      
      sub_Traj_df = rbind(sub_Traj_df, cur_traj_df)
      
    }
    
    sub_Traj_df$Sample = stringr::str_c("Sample", j)
    
    Full_Traj_df = rbind(Full_Traj_df, sub_Traj_df)
    
  }
  
  ggplot2::ggplot(data = Full_Traj_df, ggplot2::aes(x = t, y = Trajectory, fill = Sample, color = Model)) + ggplot2::geom_line() + ggplot2::scale_fill_manual(values = rep("black",length(unique(Full_Traj_df$Sample)))) + ggplot2::guides(fill = F)
  
}