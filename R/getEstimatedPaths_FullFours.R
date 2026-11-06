#' Get Estimated Path using the Full Fours Method
#'
#' @description
#' A hybrid of the ByFours and FullRegression path methods. The ByFours
#' algorithm is run first to quickly search for the spline nodes. A full cubic
#' spline regression is then fit on the selected nodes, which gives the stable
#' estimates and the estimable covariance structure of the FullRegression
#' method at close to the computational cost of ByFours.
#'
#' Because the final fit is a FullRegression fit, the result can be used
#' anywhere a FullRegression or Projection path is required, for example
#' [getEstimatedTrajectories_OneShot()], [coverageRatePath()],
#' [getRandomEstPathPosition()] or `Random_Path = TRUE`.
#'
#' @details
#' How the nodes are chosen depends on `Model_Selection_Type`:
#' * `"BIC"` or `"AIC"`: the criterion is applied to the ByFours candidates
#'   directly, and one full regression is fit on the chosen nodes. This is the
#'   fast option.
#' * `"BICVel"`, `"VarPos"`, `"VarVel"` or `"VarPosVel"`: these criteria need a
#'   covariance structure, so a full regression is fit on the nodes of every
#'   ByFours candidate and the criterion is applied to those fits. This is
#'   slower (one full regression per candidate) but still avoids the
#'   FullRegression node search.
#'
#' `patience` stops the search once that many candidates in a row have failed
#' to improve the criterion. With `"BIC"` or `"AIC"` it stops the ByFours node
#' search; with the covariance-based criteria the (cheap) node search always
#' runs in full and `patience` stops the full regressions on its candidates.
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param nsplits a numeric scalar. The maximum number of splits the ByFours step can split into.
#' @param V_Smooth a boolean. TRUE if the final spline should have C2 smoothness. FALSE if C1.
#' The ByFours node search always uses C1 smoothness.
#' @param Model_Selection_Type a string. The model criterion.
#' Must pick from the following: BIC, AIC, BICVel, VarPos, VarVel, VarPosVel. See Details.
#' @param patience NULL or a positive integer. If NULL (the default), every candidate up to `nsplits` is considered.
#' Otherwise, stop after this many consecutive candidates that do not improve the selection criterion. See Details.
#' This is faster, but it can miss a better model that only appears after a long run without improvement.
#'
#' @returns an EstimatedPath object, with the same elements as one from
#' [getEstimatedPaths_FullRegression()].
#' @export
#'
#' @examples
#' ExEstPath_FullFours <- getEstimatedPaths_FullFours(ExData, nsplits = 100)
#' coverageRatePath(ExEstPath_FullFours, ExData, CL = 0.95)
#'
getEstimatedPaths_FullFours = function(Data, nsplits, V_Smooth = F, Model_Selection_Type = "BIC", patience = NULL){

  checkPathSelectionType(Model_Selection_Type)
  checkPatience(patience)

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c('t', stringr::str_c("X",1:n_dim))

  fitFull = function(time_splits){
    tryCatch(expr = runEstPathStep_FullRegression(Data = Data, time_splits = time_splits, V_Smooth = V_Smooth),
             error = function(e){NULL})
  }

  fast_type = Model_Selection_Type %in% c("BIC", "AIC")

  # Step 1: ByFours node search (patience applies here only for BIC/AIC)
  candidates = getEstimatedPaths_ByFours(Data = Data, nsplits = nsplits, return_all = T,
                                         Model_Selection_Type = if(fast_type) Model_Selection_Type else "BIC",
                                         patience = if(fast_type) patience else NULL)

  if(fast_type){

    # Step 2: choose among the ByFours candidates, then one full regression
    time_splits = getBestSplitsPath(candidates, show = F, V_Smooth = F, type = Model_Selection_Type)$TimeSplits

    EstPath = fitFull(time_splits)

    if(is.null(EstPath)){
      stop("The full regression on the ByFours nodes (", length(time_splits) - 1,
           " regions) was singular. Try a smaller nsplits.", call. = FALSE)
    }

  } else{

    # Step 2: a full regression on each candidate, then choose among those
    tracker = pathPatienceTracker(patience, Model_Selection_Type, V_Smooth)
    full_fits = list()

    for(cand in candidates){
      fit = fitFull(cand$TimeSplits)
      if(is.null(fit)){
        next
      }
      full_fits[[length(full_fits) + 1]] = fit
      if(!is.null(patience) && tracker$add(fit)){
        break
      }
    }

    if(length(full_fits) == 0){
      stop("The full regression was singular for every ByFours candidate. Try a smaller nsplits.", call. = FALSE)
    }

    EstPath = if(is.null(patience)){
      getBestSplitsPath(full_fits, show = F, V_Smooth = V_Smooth, type = Model_Selection_Type)
    } else{
      selectBestPath(full_fits, tracker$values(), Model_Selection_Type)
    }

  }

  EstPath

}
