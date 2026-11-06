#' visualizeEstPath
#'
#' @description
#' Plot an estimated path with the data and, when the path has a covariance
#' structure (FullFours, FullRegression or Projection), its uncertainty.
#'
#' @details
#' The interval for dimension \eqn{d} at time \eqn{t} uses the variance
#' \eqn{b(t)^T \Sigma_d b(t)}, where \eqn{b(t)} is the cubic basis and
#' \eqn{\Sigma_d} the coefficient covariance. This gives a confidence interval
#' for the path itself. With `P_Int = TRUE` the observation noise
#' \eqn{\sigma^2_d} is added, which gives a prediction interval for a new
#' observation: the interval that [coverageRatePath()] checks the data against.
#'
#' How the uncertainty is drawn depends on the axes:
#' * If one axis is time (`"t"`), the interval for the other variable is drawn
#'   as a shaded band.
#' * If both axes are positions, an ellipse is drawn at `n_ellipses` evenly
#'   spaced times. The dimensions are estimated separately, so the ellipses are
#'   axis-aligned. They are joint regions for the two coordinates at level `CL`
#'   (radius \eqn{\sqrt{\chi^2_{2, CL}}} standard deviations), so they are wider
#'   than the per-dimension intervals.
#'
#' Paths without a covariance structure (ByFours) are drawn without intervals.
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#' @param EstPath an EstimatedPath object.
#' @param t_grid_size a positive numeric scalar. The time grid size of the plotting.
#' @param x_var a string. The variable name to plot on the x-axis: "t", "X1", "X2", ...
#' @param y_var a string. The variable name to plot on the y-axis.
#' @param CL a real number between 0 and 1. The confidence level of the intervals.
#' @param P_Int a boolean. TRUE for prediction intervals (path uncertainty plus observation noise),
#' FALSE for confidence intervals on the path.
#' @param n_ellipses a positive integer. The number of ellipses drawn when both axes are positions.
#'
#' @returns a ggplot2 object.
#' @export
#'
#' @examples
#' # Position plot with 95% prediction ellipses
#' visualizeEstPath(ExData, ExEstPath_FullFours, x_var = 'X1', y_var = 'X2')
#'
#' # One dimension over time with a 95% confidence band for the path
#' visualizeEstPath(ExData, ExEstPath_FullFours, x_var = 't', y_var = 'X1', P_Int = FALSE)
#'
visualizeEstPath = function(Data, EstPath, t_grid_size = 0.01, x_var = "t", y_var = "X1", CL = 0.95, P_Int = T, n_ellipses = 20){

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  dim_names = stringr::str_c("X", 1:n_dim)
  names(Data) = c("t", dim_names)

  if(!all(c(x_var, y_var) %in% names(Data)) || x_var == y_var){
    stop("x_var and y_var must be two different names among: ", paste(names(Data), collapse = ", "), call. = FALSE)
  }

  if(length(CL) != 1 || !is.numeric(CL) || CL <= 0 || CL >= 1){
    stop("CL must be a number between 0 and 1.", call. = FALSE)
  }

  TimeSplits = EstPath$TimeSplits
  has_cov = !is.null(EstPath$Cov) && !is.null(EstPath$Sigma2)

  # Path estimate and per-dimension standard deviation on a time grid
  pathWithSD = function(t){
    path = evaluatePathAt(t, EstPath)
    sd = matrix(NA_real_, length(t), n_dim)
    if(has_cov){
      B = cbind(path$s^3, path$s^2, path$s, 1)
      for(r in unique(path$idx)){
        sel = path$idx == r
        for(d in seq_len(n_dim)){
          Bs = B[sel, , drop = FALSE]
          v = rowSums((Bs %*% EstPath$Cov[[r]][[d]]) * Bs) + P_Int*EstPath$Sigma2[[r]][d]
          sd[sel, d] = sqrt(pmax(v, 0))
        }
      }
    }
    colnames(path$pos) = dim_names
    colnames(sd) = dim_names
    list(pos = path$pos, sd = sd)
  }

  t_seq = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)
  grid = pathWithSD(t_seq)
  Full_Path_df = data.frame(t = t_seq, grid$pos)

  x_sym = dplyr::sym(x_var)
  y_sym = dplyr::sym(y_var)

  plot = ggplot2::ggplot() +
    ggplot2::geom_point(data = Data, ggplot2::aes(x = !!x_sym, y = !!y_sym), alpha = 0.15)

  interval_label = stringr::str_c(round(100*CL), "% ", if(P_Int) "prediction" else "confidence")

  if(has_cov){

    if(x_var == "t" || y_var == "t"){

      # Shaded band for the non-time variable
      dim_var = if(x_var == "t") y_var else x_var
      z = stats::qnorm((1 + CL)/2)
      Full_Path_df$lower = grid$pos[, dim_var] - z*grid$sd[, dim_var]
      Full_Path_df$upper = grid$pos[, dim_var] + z*grid$sd[, dim_var]

      if(x_var == "t"){
        plot = plot + ggplot2::geom_ribbon(data = Full_Path_df, ggplot2::aes(x = t, ymin = lower, ymax = upper), fill = "red", alpha = 0.3)
      } else{
        plot = plot + ggplot2::geom_ribbon(data = Full_Path_df, ggplot2::aes(y = t, xmin = lower, xmax = upper), fill = "red", alpha = 0.3, orientation = "y")
      }

      plot = plot + ggplot2::labs(caption = stringr::str_c("Shaded: ", interval_label, " interval"))

    } else{

      # Axis-aligned joint ellipses at evenly spaced times
      t_ell = seq(min(TimeSplits), max(TimeSplits), length.out = n_ellipses)
      ell = pathWithSD(t_ell)
      radius = sqrt(stats::qchisq(CL, df = 2))
      angle = seq(0, 2*pi, length.out = 61)

      Ellipse_df = do.call(rbind, lapply(seq_along(t_ell), function(i){
        data.frame(id = i,
                   x = ell$pos[i, x_var] + radius*ell$sd[i, x_var]*cos(angle),
                   y = ell$pos[i, y_var] + radius*ell$sd[i, y_var]*sin(angle))
      }))

      plot = plot + ggplot2::geom_polygon(data = Ellipse_df, ggplot2::aes(x = x, y = y, group = id), fill = "red", alpha = 0.15, colour = "red", linewidth = 0.3) +
        ggplot2::labs(caption = stringr::str_c("Ellipses: joint ", interval_label, " regions"))

    }
  }

  plot + ggplot2::geom_path(data = Full_Path_df, ggplot2::aes(x = !!x_sym, y = !!y_sym))

}
