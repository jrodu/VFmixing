#' visualizeEstPath
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#' @param EstPath an EstimatedPath object.
#' @param t_grid_size a positive numeric scalar. The time grid size of the plotting.
#' @param x_var a string. The variable name to plot on the x-axis.
#' @param y_var a string. The variable name to plot on the y-axis.
#'
#' @returns a ggplot2 object.
#' @export
#' 
#' @examples
#' visualizeEstPath(ExData, ExEstPath_ByFours, t_grid_size = 0.01, 
#'                  x_var = 'X1',y_var = 'X2')
#'       
visualizeEstPath = function(Data, EstPath, t_grid_size = 0.01, x_var = "t", y_var = "X1"){
  Data = data.frame(Data)
  names(Data) = c("t",stringr::str_c("X",1:(ncol(Data)-1)))
  TimeSplits = EstPath$TimeSplits
  
  t_seq = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)
  err_dim = as.numeric(substr(y_var, 2, nchar(y_var)))
  
  evaluateIndividualPath = function(t){
    cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
    cur_pos = evaluateCubic(t, EstPath$Path[[cur_path_index]], startTime = TimeSplits[cur_path_index])
    cur_difft = t - TimeSplits[cur_path_index]
    trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)
    
    upper = cur_pos[err_dim]
    lower = cur_pos[err_dim]
  
    c(t, cur_pos, upper, lower)
  }
  
  Full_Path_df = data.frame(t(apply(matrix(t_seq), 1, evaluateIndividualPath)))
  
  names(Full_Path_df) = c('t',stringr::str_c("X",1:(ncol(Data)-1)), 'upper','lower')
  
  ggplot2::ggplot(data = Full_Path_df, ggplot2::aes(x = !!dplyr::sym(x_var), y = !!dplyr::sym(y_var))) + ggplot2::geom_path() + ggplot2::geom_point(data = Data, ggplot2::aes(x = !!dplyr::sym(x_var), y = !!dplyr::sym(y_var)), alpha = 0.05) + 
    ggplot2::geom_ribbon(ggplot2::aes(x=!!dplyr::sym(x_var), y=!!dplyr::sym(y_var), ymax=upper, ymin=lower), alpha=0.5, fill = 'red')
  
}
