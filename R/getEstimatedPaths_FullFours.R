#' Get Estimated Path using the Full Fours Method
#'
#' @description
#' A hybrid of the ByFours and FullRegression path methods. The ByFours
#' algorithm is run first to quickly choose the spline nodes (selected by BIC).
#' A single full cubic spline regression is then fit on those nodes, which
#' gives the stable estimates and the estimable covariance structure of the
#' FullRegression method at close to the computational cost of ByFours.
#'
#' Because the final fit is a FullRegression fit, the result can be used
#' anywhere a FullRegression or Projection path is required, for example
#' [getEstimatedTrajectories_OneShot()], [coverageRatePath()],
#' [getRandomEstPathPosition()] or `Random_Path = TRUE`.
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param nsplits a numeric scalar. The maximum number of splits the ByFours step can split into.
#' @param V_Smooth a boolean. TRUE if the final spline should have C2 smoothness. FALSE if C1.
#' The ByFours node search always uses C1 smoothness.
#'
#' @returns an EstimatedPath object, with the same elements as one from
#' [getEstimatedPaths_FullRegression()].
#' @export
#'
#' @examples
#' ExEstPath_FullFours <- getEstimatedPaths_FullFours(ExData, nsplits = 100)
#' coverageRatePath(ExEstPath_FullFours, ExData, CL = 0.95)
#'
getEstimatedPaths_FullFours = function(Data, nsplits, V_Smooth = F){

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))

  # Step 1: ByFours node search
  ByFours_Path = getEstimatedPaths_ByFours(Data = Data, nsplits = nsplits)

  time_splits = ByFours_Path$TimeSplits

  # Step 2: one full regression on the selected nodes
  EstPath = tryCatch(expr = runEstPathStep_FullRegression(Data = Data, time_splits = time_splits, V_Smooth = V_Smooth),
                     error = function(e){NULL})

  if(is.null(EstPath)){
    stop("The full regression on the ByFours nodes (", length(time_splits) - 1,
         " regions) was singular. Try a smaller nsplits.", call. = FALSE)
  }

  EstPath

}
