#' getPropagationList
#'
#' @param Path an EstimatedPath object.
#' @param TrajList a list of EstimatedTrajectory objects.
#' @param baseVectorFields a baseVectorFields function.
#' @param prop_t_grid_size the time grid size of the propagation points.
#' @param n_prop_steps the number of intermediate RungeKutta steps for each propagation.
#'
#' @returns a PropagationList
#' @export
#'
#' @examples
#' getPropagationList(Path = ExEstPath_ByFours, 
#'              TrajList = list(ExEstTraj_ByFours, ExEstTraj_Optimization), 
#'              baseVectorFields = baseVectorFields, 
#'             prop_t_grid_size = 0.1, n_prop_steps = 10)
#' 
getPropagationList = function(Path, TrajList, baseVectorFields, prop_t_grid_size, n_prop_steps){
  
  PropagationList = list()
  n_methods = length(TrajList)
  
  t_seq = seq(min(Path$TimeSplits), max(Path$TimeSplits), by = prop_t_grid_size)
  
  PathData = data.frame(cbind(t_seq, getEstPathPosition(t_seq, Path)))
  names(PathData) = c("Time","Longitude","Latitude")
  
  RK_ints = diff(PathData$Time)
  
  # 1. Initialize a list to store the final estimated positions for each method.
  est_pos_list = replicate(n_methods, list())
  
  PropagationList[[2]] = list()
  
  for(k in 1:length(RK_ints)){
    
    # Create temporary lists to hold results for the current step 'k'.
    cur_RK_list = list()
    
    # 2. Loop through each trajectory method.
    for (j in 1:n_methods) {
      cur_RK = RungeKutta(startTime = PathData$Time[k], startPos = c(PathData[k,2],PathData[k,3]),
                          baseVectorFields = baseVectorFields, 
                          TrajList = TrajList[[j]]$Traj, 
                          TrajTimeSplits = TrajList[[j]]$TimeSplits, 
                          endTime = PathData$Time[k+1], 
                          t_step = RK_ints[k]/n_prop_steps)
      
      # Add the result to our list for this step
      cur_RK_list[[j]] = cur_RK
      
      # Append the final position from this RK run to the corresponding list
      est_pos_list[[j]][[k]] = cur_RK[nrow(cur_RK),]
    }
    
    # 3. Combine the results for this k-step dynamically.
    # Start with the first trajectory's results (which includes the 'time' column).
    combined_RK_step = cur_RK_list[[1]]
    # If there's more than one method, cbind the remaining results without their time columns.
    if (n_methods > 1) {
      for (j in 2:n_methods) {
        combined_RK_step = cbind(combined_RK_step, cur_RK_list[[j]][, -1, drop = FALSE])
      }
    }
    
    PropagationList[[2]][[k]] = combined_RK_step
    
    # 4. Create column names dynamically.
    names(PropagationList[[2]][[k]]) = c('time', paste0(c("LonM", "LatM"), rep(1:n_methods, each = 2)))
    
    svMisc::progress(k, length(RK_ints))
  }
  
  # 5. Combine the lists of final positions into data frames.
  est_pos_dfs = lapply(est_pos_list, function(x) do.call(rbind, x))
  
  # 6. Dynamically calculate the final differences.
  # Calculate the difference for each method compared to the true path.
  diff_list = lapply(est_pos_dfs, function(df) {
    PathData[-1, -1] - df[, -1] 
  })
  
  # Combine the 'Time' column from the first method's estimates with all differences.
  PropagationList[[1]] = cbind(est_pos_dfs[[1]][, 1, drop = FALSE], do.call(cbind, diff_list))
  
  # 7. Set the final column names dynamically.
  colnames(PropagationList[[1]]) = c('Time', paste0(c("LonM", "LatM"), rep(1:n_methods, each=2)))
  
  PropagationList
}
