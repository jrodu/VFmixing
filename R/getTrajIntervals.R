#' Gridded Trajectory Estimates and Confidence Intervals
#'
#' @description
#' Evaluates the estimated trajectories and their pointwise confidence
#' intervals on a time grid, one data frame per trajectory: the values
#' [visualizeEstTraj()] plots and [coverageRateTraj()] checks.
#'
#' @details
#' At each time \eqn{t} the interval is
#' \eqn{\hat\beta(t) \pm z \sqrt{s_w^2(t) + V_{sel}(t)}}, where \eqn{s_w^2(t)} is
#' the within-node variance from the trajectory's covariance (`Cov`) and
#' \eqn{V_{sel}(t)} the node-selection variance from its bootstrap conditional
#' means (`SelectionMeans`, from [getEstimatedTrajectories_Full()] with
#' `bootstrap = TRUE`; zero without them). A trajectory with no covariance
#' (for example from [getEstimatedTrajectories_ByFours()]) gets the estimates
#' with `NA` bounds.
#'
#' @param EstTraj an EstimatedTrajectory object.
#' @param t_grid_size a positive number. The spacing of the time grid, in the data's time units.
#' @param CL a number between 0 and 1. The confidence level of the intervals.
#' @param labels NULL or a character vector with one name per trajectory, used to name the list. NULL gives
#' "Trajectory1", "Trajectory2", ...
#'
#' @returns a named list with one data frame per trajectory, with columns `t`, `estimate`, `lower` and `upper`.
#' @export
#'
#' @examples
#' intervals = getTrajIntervals(ExEstTraj_OneShot, t_grid_size = 0.1, labels = c("Field A", "Field B"))
#' head(intervals[["Field A"]])
#'
getTrajIntervals = function(EstTraj, t_grid_size = 0.01, CL = 0.95, labels = NULL){

  checkIntervalArgs(CL, t_grid_size)

  n_models = nrow(EstTraj$Traj[[1]])
  labels = trajLabels(labels, n_models)

  t = seq(min(EstTraj$TimeSplits), max(EstTraj$TimeSplits), by = t_grid_size)
  estimate = evaluateTrajAt(t, EstTraj)

  if(is.null(EstTraj$Cov)){
    half_width = matrix(NA_real_, length(t), n_models)
  } else{
    variance = trajWithinVariance(EstTraj, t) + trajSelectionVariance(EstTraj, t)
    half_width = stats::qnorm((1 + CL)/2)*sqrt(pmax(variance, 0))
  }

  out = lapply(seq_len(n_models), function(m){
    data.frame(t = t, estimate = estimate[, m], lower = estimate[, m] - half_width[, m], upper = estimate[, m] + half_width[, m])
  })
  names(out) = labels
  out

}

# Trajectory labels: the given ones (one per trajectory), or "Trajectory1", ...
trajLabels = function(labels, n_models){

  if(is.null(labels)){
    return(stringr::str_c("Trajectory", seq_len(n_models)))
  }
  if(length(labels) != n_models || anyDuplicated(labels) || any(is.na(labels))){
    stop("labels must give ", n_models, " distinct names, one per trajectory.", call. = FALSE)
  }
  as.character(labels)

}
