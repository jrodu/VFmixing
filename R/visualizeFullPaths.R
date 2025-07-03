#' Visualize Path Functions
#'
#' @description
#' A function to visualize a series of path functions. 
#' 
#'
#' @param PathList a list of lists. Each element of the outer list corresponds to a list of cubic matrix coefficients associated with individual cubic functions.
#' @param TimeSplitsList a list of numeric vectors. Each element corresponds to nodes of the cubic spline in the current element of PathList.
#' @param t_grid_size a non-negative numeric scalar. The time grid size of the plotting.
#'
#' @returns a ggplot2 object. 
#' @export
#'
#' @examples
#' 
#' Path <- list(matrix(c(0,0,0,1,
#'                       0,0,0,1), nrow = 2, byrow = TRUE))
#' TimeSplits <- c(0,1)
#' 
#' visualizeFullPaths(PathList = list(Path), TimeSplitsList = list(TimeSplits), t_grid_size = 0.01)
visualizeFullPaths = function(PathList, TimeSplitsList, t_grid_size){
  Full_Path_df = data.frame()
  for(j in 1:length(PathList)){
    
    cur_PathList = PathList[[j]]
    TimeSplits = TimeSplitsList[[j]]
    
    sub_Path_df = data.frame()
    
    for(i in 1:length(cur_PathList)){
      cur_t_seq = seq(TimeSplits[i], TimeSplits[i+1], by = t_grid_size) 
      
      cur_path_df = data.frame(cbind(cur_t_seq,
                                     evaluateCubic(t = cur_t_seq, Cubic = cur_PathList[[i]], startTime = TimeSplits[i])))
      names(cur_path_df) = c("t",stringr::str_c("X",1:(ncol(cur_path_df)-1)))
      cur_path_df = data.frame(apply(cur_path_df, 2, function(x) as.numeric(as.character(x))))
      sub_Path_df = rbind(sub_Path_df, cur_path_df)
      
    }
    
    sub_Path_df$Sample = stringr::str_c("Sample", j)
    
    Full_Path_df = rbind(Full_Path_df, sub_Path_df)
    
  }
  
  ggplot2::ggplot(data = Full_Path_df, ggplot2::aes(x = X1, y = X2, fill = Sample)) + ggplot2::geom_path(ggplot2::aes(text = paste0("Time: ",t,"\n Longitude: ",X1,"\n Latitude: ",X2))) + ggplot2::scale_fill_manual(values = rep("black",length(unique(Full_Path_df$Sample)))) + ggplot2::guides(fill = F)
  
}