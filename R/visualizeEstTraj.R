#' visualizeEstTraj
#'
#' @param EstTraj an EstimatedTrajectory object. 
#' @param TrueTraj a list of numeric matrices. The list of coefficient matrices of a true trajectory. Not required.
#' @param TrueTrajTimeSplits a numeric vector. The nodes of the true trajectory in addition with the starting and ending time. Not required.
#' @param EstPath an EstimatedPath object. The path used to estimate the trajectories.
#' @param baseVectorFields a baseVectorFieldFunction.
#' @param t_grid_size a positive numeric scalar. The time grid size of the plotting.
#' @param CL a real number between 0 and 1. The confidence level of the intervals to plot. Only used if EstTraj was estimated using the Optimization method.
#' @param P_Int a boolean. TRUE if the prediction interval should be plotted. FALSE if the confidence interval on the mean should be plotted.
#' Only used if EstTraj was estimated using the Optimization method.
#' @param Cos_Cutoff a real number between 0 and 1. The cosine value to determine the cutoff for when vector fields are too colinear. 
#'
#' @returns a ggplot2 object.
#' @export
#' 
#' @examples
#' # visualizeEstTraj(EstTraj = ExEstTraj_Optimization, TrueTraj = ExTraj$Traj, 
#' #                  TrueTrajTimeSplits = ExTraj$TimeSplits, 
#' #                  EstPath = ExEstPath_Projection, 
#' #                  baseVectorFields = baseVectorFields, t_grid_size = 0.01, 
#' #                  Cos_Cutoff = 0.91)
#'   
visualizeEstTraj = function(EstTraj, TrueTraj, TrueTrajTimeSplits, EstPath, baseVectorFields, t_grid_size = 0.01, CL = 0.95, P_Int = T, Cos_Cutoff = 0.91){
  
  n_models = nrow(EstTraj$Traj[[1]])
  n_dim = nrow(EstPath$Path[[1]])
  TimeSplits_Traj = EstTraj$TimeSplits
  TimeSplits_Path = EstPath$TimeSplits
  Cov = EstTraj$Cov
  
  t_seq = seq(min(TimeSplits_Traj), max(TimeSplits_Traj), by = t_grid_size)
  
  evaluateIndividualEstTraj = function(t){
    cur_traj_index = min(c(max(c(sum((TimeSplits_Traj - t)<0),1)), length(TimeSplits_Traj)-1))
    cur_path_index = min(c(max(c(sum((TimeSplits_Path - t)<0),1)), length(TimeSplits_Path)-1))
    cur_traj = evaluateCubic(t, EstTraj$Traj[[cur_traj_index]], startTime = TimeSplits_Traj[cur_traj_index])
    cur_pos = evaluateCubic(t, EstPath$Path[[cur_path_index]], TimeSplits_Path[cur_path_index])
    BaseVel = baseVectorFields(t, cur_pos)
    
    if(n_models > 1){
      model_comb = t(utils::combn(1:n_models, m = 2))
      
      all_angles = c()
      
      for(i in 1:nrow(model_comb)){
        
        all_angles = c(all_angles, abs(sum(BaseVel[,model_comb[i,1]] * BaseVel[,model_comb[i,2]]) / (sqrt(sum(BaseVel[,model_comb[i,1]]^2))*sqrt(sum(BaseVel[,model_comb[i,2]]^2)))))
        
      }
    } else{
      all_angles = c(0)
    }
    
    
    if(is.null(Cov)){
      upper = cur_traj
      lower = cur_traj
    } else{
      cur_difft = t - TimeSplits_Traj[cur_traj_index]
      trans_mat = matrix(rep(c(cur_difft^3,cur_difft^2,cur_difft,1), n_models), nrow = 4)
      
      cur_y_var = apply(matrix(1:n_models), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_traj_index]][1:4 + (x-1)*4, 1:4 + (x-1)*4] %*% trans_mat)[1,1] + (EstTraj$Sigma2[[cur_traj_index]])*P_Int})
      
      upper = cur_traj + stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
      lower = cur_traj - stats::qnorm((1+CL)/2)*sqrt(cur_y_var)
    }
    
    cbind(t(t(rep(t, n_models))), t(cur_traj), t(upper), t(lower), t(t(stringr::str_c("Trajectory", 1:n_models))), t(t(rep(max(all_angles) < Cos_Cutoff, n_models))))
  }
  
  evaluateIndividualTrueTraj = function(t){
    
    cur_true_traj_index = min(c(max(c(sum((TrueTrajTimeSplits - t)<0),1)), length(TrueTrajTimeSplits)-1))
    cur_traj = evaluateCubic(t, TrueTraj[[cur_true_traj_index]], startTime = TrueTrajTimeSplits[cur_true_traj_index])
    
    cbind(t(t(rep(t, n_models))), t(cur_traj), t(t(stringr::str_c("Trajectory", 1:n_models))))
    
  }
  
  Full_EstTraj_df = data.frame(do.call("rbind", lapply(as.list(t_seq), evaluateIndividualEstTraj)))
  names(Full_EstTraj_df) = c('t','Trajectory', 'upper','lower','Model','Angle')
  
  if(!missing(TrueTraj)){
    Full_TrueTraj_df = data.frame(do.call("rbind", lapply(as.list(t_seq), evaluateIndividualTrueTraj)))
    names(Full_TrueTraj_df) = c('t','Trajectory','Model')
    Full_TrueTraj_df$t = as.numeric(Full_TrueTraj_df$t)
    Full_TrueTraj_df$Trajectory = as.numeric(Full_TrueTraj_df$Trajectory)
  }

  Full_EstTraj_df$t = as.numeric(Full_EstTraj_df$t)
  Full_EstTraj_df$Trajectory = as.numeric(Full_EstTraj_df$Trajectory)
  Full_EstTraj_df$upper = as.numeric(Full_EstTraj_df$upper)
  Full_EstTraj_df$lower = as.numeric(Full_EstTraj_df$lower)
  Full_EstTraj_df$Angle[Full_EstTraj_df$Angle == "TRUE"] = stringr::str_c("Cos(", paste("\U03B8"),") < ", round(Cos_Cutoff,2))
  Full_EstTraj_df$Angle[Full_EstTraj_df$Angle == "FALSE"] = stringr::str_c("Cos(", paste("\U03B8"),") > ", round(Cos_Cutoff,2))
  
  n_colors <- length(unique(Full_EstTraj_df$Model))
  trajectory_colors <- stats::setNames(RColorBrewer::brewer.pal(max(3, n_colors), "Set1")[1:n_colors], 
                                unique(Full_EstTraj_df$Model))
  
  angle_vT = stringr::str_c("Cos(", paste("\U03B8"),") < ", round(Cos_Cutoff,2))
  angle_vF = stringr::str_c("Cos(", paste("\U03B8"),") > ", round(Cos_Cutoff,2))
  
  background_colors <- c("green","red")
  names(background_colors) = c(angle_vT, angle_vF)
  
  
 plot = ggplot2::ggplot() + 
   ggplot2::geom_path(data = Full_EstTraj_df, ggplot2::aes(x = t, y = Trajectory, color = Model), linetype = 'solid') + 
   ggplot2::geom_rect(data = Full_EstTraj_df, ggplot2::aes(xmin = t, xmax = dplyr::lead(t, default = max(t)), ymin = -Inf, ymax = Inf, fill = factor(Angle)), alpha = 0.1) +  
   ggplot2::geom_ribbon(data = Full_EstTraj_df, ggplot2::aes(x = t, ymin = lower, ymax = upper, fill = Model), alpha = 0.5) +
   ggplot2::scale_color_manual(values = trajectory_colors) +
   ggplot2::scale_fill_manual(
      name = "Model Confidence \n & Model Angle",
      values = c(trajectory_colors, background_colors), 
      guide = ggplot2::guide_legend(order = 1)
    ) +
   ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1), 
      color = ggplot2::guide_legend(order = 2)
    )
 
 if(!missing(TrueTraj)){
   plot = plot + ggplot2::geom_path(data = Full_TrueTraj_df, ggplot2::aes(x = t, y = Trajectory, color = Model), linetype = 'dashed')
 }
 
 plot
  
}
