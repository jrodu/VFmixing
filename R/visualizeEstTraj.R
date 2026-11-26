#' visualizeEstTraj
#'
#' @description
#' Plots the estimated trajectories with their pointwise confidence bands (from
#' [getTrajIntervals()]), optionally with the true trajectories and with the
#' times shaded where the base vector fields are close to collinear.
#'
#' @param EstTraj an EstimatedTrajectory object.
#' @param TrueTraj a list of numeric matrices. The list of coefficient matrices of a true trajectory. Not required.
#' @param TrueTrajTimeSplits a numeric vector. The nodes of the true trajectory in addition with the starting and ending
#' time. Needed with TrueTraj.
#' @param EstPath an EstimatedPath object, the path used to estimate the trajectories. Not required: with
#' `baseVectorFields`, it is used to shade the times where the fields along the path are close to collinear.
#' @param baseVectorFields a baseVectorFieldFunction. Not required: see `EstPath`.
#' @param t_grid_size a positive numeric scalar. The time grid size of the plotting.
#' @param CL a real number between 0 and 1. The confidence level of the intervals to plot. Only used if EstTraj has a
#' covariance structure (OneShot). The bands are confidence intervals from the trajectory coefficient covariance,
#' plus the node-selection variance if EstTraj comes from [getEstimatedTrajectories_Full()] with the bootstrap. There
#' is no prediction interval: the trajectories are never observed directly.
#' @param Cos_Cutoff a real number between 0 and 1. The cosine value to determine the cutoff for when vector fields
#' are too colinear. Only used with `EstPath` and `baseVectorFields`.
#' @param labels NULL or a character vector with one name per trajectory, shown in the legend. NULL gives
#' "Trajectory1", "Trajectory2", ...
#'
#' @returns a ggplot2 object.
#' @export
#'
#' @examples
#' visualizeEstTraj(EstTraj = ExEstTraj_OneShot, t_grid_size = 0.05, labels = c("Field A", "Field B"))
#'
#' # With the truth and the collinearity shading
#' # visualizeEstTraj(EstTraj = ExEstTraj_OneShot, TrueTraj = ExTraj$Traj,
#' #                  TrueTrajTimeSplits = ExTraj$TimeSplits,
#' #                  EstPath = ExEstPath_FullFours,
#' #                  baseVectorFields = baseVectorFields, t_grid_size = 0.01,
#' #                  Cos_Cutoff = 0.91)
#'
visualizeEstTraj = function(EstTraj, TrueTraj = NULL, TrueTrajTimeSplits = NULL, EstPath = NULL, baseVectorFields = NULL,
                            t_grid_size = 0.01, CL = 0.95, Cos_Cutoff = 0.91, labels = NULL){

  n_models = nrow(EstTraj$Traj[[1]])
  labels = trajLabels(labels, n_models)

  if(!is.null(TrueTraj) && is.null(TrueTrajTimeSplits)){
    stop("TrueTrajTimeSplits is needed with TrueTraj.", call. = FALSE)
  }

  shade_angles = !is.null(EstPath) && !is.null(baseVectorFields)

  # The estimates and bands (the same values getTrajIntervals() returns)
  intervals = getTrajIntervals(EstTraj, t_grid_size = t_grid_size, CL = CL, labels = labels)
  t_seq = intervals[[1]]$t

  Full_EstTraj_df = do.call(rbind, lapply(labels, function(label){
    d = intervals[[label]]
    # No covariance: a zero-width band at the estimate
    no_band = is.na(d$lower)
    d$lower[no_band] = d$estimate[no_band]
    d$upper[no_band] = d$estimate[no_band]
    data.frame(t = d$t, Trajectory = d$estimate, upper = d$upper, lower = d$lower, Model = label)
  }))
  Full_EstTraj_df$Model = factor(Full_EstTraj_df$Model, levels = labels)

  n_colors = length(labels)
  trajectory_colors = stats::setNames(RColorBrewer::brewer.pal(max(3, n_colors), "Set1")[1:n_colors], labels)

  plot = ggplot2::ggplot() +
    ggplot2::geom_path(data = Full_EstTraj_df, ggplot2::aes(x = t, y = Trajectory, color = Model), linetype = 'solid')

  fill_values = trajectory_colors
  fill_name = "Model Confidence"

  if(shade_angles){

    # The largest |cos(angle)| between any two fields along the estimated path
    positions = getEstPathPosition(t_seq, EstPath)
    max_cos = vapply(seq_along(t_seq), function(i){
      BaseVel = baseVectorFields(t_seq[i], positions[i, ])
      if(n_models < 2){
        return(0)
      }
      pairs = utils::combn(n_models, 2)
      max(apply(pairs, 2, function(p){
        abs(sum(BaseVel[, p[1]]*BaseVel[, p[2]]))/(sqrt(sum(BaseVel[, p[1]]^2))*sqrt(sum(BaseVel[, p[2]]^2)))
      }))
    }, numeric(1))

    angle_vT = stringr::str_c("Cos(", paste("\U03B8"), ") < ", round(Cos_Cutoff, 2))
    angle_vF = stringr::str_c("Cos(", paste("\U03B8"), ") > ", round(Cos_Cutoff, 2))
    angle_df = data.frame(t = t_seq, Angle = ifelse(max_cos < Cos_Cutoff, angle_vT, angle_vF))

    background_colors = c("green", "red")
    names(background_colors) = c(angle_vT, angle_vF)
    fill_values = c(trajectory_colors, background_colors)
    fill_name = "Model Confidence \n & Model Angle"

    plot = plot +
      ggplot2::geom_rect(data = angle_df, ggplot2::aes(xmin = t, xmax = dplyr::lead(t, default = max(t)), ymin = -Inf, ymax = Inf,
                                                       fill = factor(Angle)), alpha = 0.1)
  }

  plot = plot +
    ggplot2::geom_ribbon(data = Full_EstTraj_df, ggplot2::aes(x = t, ymin = lower, ymax = upper, fill = Model), alpha = 0.5) +
    ggplot2::scale_color_manual(values = trajectory_colors) +
    ggplot2::scale_fill_manual(name = fill_name, values = fill_values, guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = ggplot2::guide_legend(order = 2))

  if(!is.null(TrueTraj)){
    true_values = evaluateTrajAt(t_seq, list(Traj = TrueTraj, TimeSplits = TrueTrajTimeSplits))
    Full_TrueTraj_df = do.call(rbind, lapply(seq_len(n_models), function(m){
      data.frame(t = t_seq, Trajectory = true_values[, m], Model = labels[m])
    }))
    Full_TrueTraj_df$Model = factor(Full_TrueTraj_df$Model, levels = labels)
    plot = plot + ggplot2::geom_path(data = Full_TrueTraj_df, ggplot2::aes(x = t, y = Trajectory, color = Model), linetype = 'dashed')
  }

  plot

}
