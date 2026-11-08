#' Path Model Selection
#'
#' @param EstPath a list of EstimatedPath Objects. These should all be related using the same data with an increasing number of splits.
#' @param show a boolean. TRUE if the criterion values should be displayed.
#' @param V_Smooth a boolean. TRUE if the path cubic spline had C2 smoothness. FALSE if C1.
#' @param type a string. The model criterion.
#' Must pick from the following: BIC, BICVel, AIC, VarPos, VarVel, VarPosVel.
#' BIC and AIC can be used with any path method. BICVel, VarPos, VarVel, VarPosVel need a covariance structure,
#' so they can only be used with FullRegression, Projection and FullFours.
#' @param undersmooth a boolean. If TRUE, return the first candidate with at least twice as many regions
#' as the one the criterion selects (or the last candidate, if none has that many), instead of the selected one.
#' See [getEstimatedPaths_FullFours()] for why.
#' @param best_index NULL or an integer. The index of the candidate the criterion selected, if already known
#' (for example from a search that stopped early). If NULL, it is found from the criterion.
#'
#' @returns an EstimatedPath object.
#' @export
#'
getBestSplitsPath = function(EstPath, show = F, V_Smooth = T, type = "BIC", undersmooth = F, best_index = NULL){

  checkPathSelectionType(type, EstPath[[1]])

  values = pathSelectionValues(EstPath, type, V_Smooth)

  selectBestPath(EstPath, values, type, show, undersmooth = undersmooth, best_index = best_index)

}

# Validate a path model-selection type for a given candidate path.
checkPathSelectionType = function(type, example_path = NULL){

  valid_types = c("BIC", "BICVel", "AIC", "VarPos", "VarVel", "VarPosVel")

  if(length(type) != 1 || !(type %in% valid_types)){
    stop("type must be one of: ", paste(valid_types, collapse = ", "), call. = FALSE)
  }

  if(!is.null(example_path) && type %in% c("BICVel", "VarPos", "VarVel", "VarPosVel") && is.null(example_path$Cov)){
    stop("type = \"", type, "\" needs a path with a covariance structure ",
         "(FullRegression, Projection or FullFours). Use \"BIC\" or \"AIC\" instead.", call. = FALSE)
  }

  invisible(type)
}

# The model-selection values for each candidate path. A numeric vector (smaller
# is better), except for BICVel, which gives a 2 x K matrix: the BIC in the first
# row and the 95% velocity variance in the second.
pathSelectionValues = function(EstPath, type, V_Smooth){

  if(type == "BIC"){

    getBIC = function(cur_path){

      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)
      n_param = n_dim*splineFreeParams(length(cur_path$Path), V_Smooth)

      return((n_param)*log(n_points) + n_points*log(stats::median(rowSums(cur_path$Error[,-1]))))

    }

    return(sapply(EstPath, getBIC))

  }

  if(type == "BICVel"){

    getBICVel = function(cur_path){

      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      Error = cur_path$Error
      n_param = n_dim*splineFreeParams(length(cur_path$Path), V_Smooth)

      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)

        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))

      }

      cur_VelVar = unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))

      cur_BIC = (n_param)*log(n_points) + n_points*log(stats::median(unlist(as.vector(cur_path$Error[,-1]))))

      c(cur_BIC, cur_VelVar)

    }

    return(matrix(sapply(EstPath, getBICVel), nrow = 2))

  }

  if(type == "AIC"){

    getAIC = function(cur_path){

      n_dim = nrow(cur_path$Path[[1]])
      n_points = nrow(cur_path$Error)*(ncol(cur_path$Error)-1)

      return(2*(2*n_dim*(length(cur_path$Path)+1)) + n_points*log(stats::median(unlist(as.vector(cur_path$Error[-1])))))
    }

    return(sapply(EstPath, getAIC))

  }

  if(type == "VarPos"){

    getMaxEstPosVar = function(cur_path){

      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2
      Error = cur_path$Error

      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1), nrow = 4)

        (mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]})))

      }

      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))

    }

    return(sapply(EstPath, FUN = getMaxEstPosVar, simplify = T))

  }

  if(type == "VarVel"){

    getMaxEstVelVar = function(cur_path){

      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2

      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)

        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))

      }

      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))

    }

    return(sapply(EstPath, FUN = getMaxEstVelVar, simplify = T))

  }

  if(type == "VarPosVel"){

    getMaxEstPosVelVar = function(cur_path){

      t_seq = cur_path$Error[,1]
      n_dim = nrow(cur_path$Path[[1]])
      TimeSplits = cur_path$TimeSplits
      Cov = cur_path$Cov
      Sigma2 = cur_path$Sigma2

      evaluateIndividualPath = function(t){
        cur_path_index = min(c(max(c(sum((TimeSplits - t)<0),1)), length(TimeSplits)-1))
        cur_difft = t - TimeSplits[cur_path_index]
        trans_mat = matrix(c(cur_difft^3,cur_difft^2,cur_difft,1) + c(3*cur_difft^2,2*cur_difft,1,0), nrow = 4)

        mean(apply(matrix(1:n_dim), MARGIN = 1, FUN = function(x){(t(trans_mat) %*% Cov[[cur_path_index]][[x]] %*% trans_mat)[1,1] + Sigma2[[cur_path_index]][x]}))

      }

      unname(stats::quantile(apply(matrix(t_seq), MARGIN = 1, FUN = evaluateIndividualPath), probs = c(0.95)))

    }

    return(sapply(EstPath, FUN = getMaxEstPosVelVar, simplify = T))

  }

}

# Index of the best candidate given its selection values.
pathBestIndex = function(values, type){

  if(type == "BICVel"){
    # Only consider models before the velocity variance first jumps by more
    # than a factor of e. If it never jumps, consider every model.
    vel_cutoff = which(diff(log(values[2,])) > 1)[1]

    if(is.na(vel_cutoff)){
      vel_cutoff = ncol(values)
    }

    return(which.min(values[1, 1:vel_cutoff]))
  }

  which.min(values)
}

# Pick, report and return the best candidate path. With undersmooth, return the
# first candidate with at least twice the selected number of regions instead.
selectBestPath = function(EstPath, values, type, show = F, undersmooth = F, best_index = NULL){

  if(is.null(best_index)){
    best_index = pathBestIndex(values, type)
  }

  if(undersmooth){
    best_index = undersmoothIndex(sapply(EstPath, function(p){length(p$Path)}), best_index)
  }

  reportBestRegions(length(EstPath[[best_index]]$Path))
  if(show){
    print(values)
  }

  EstPath[[best_index]]

}

# Validate a patience argument: NULL, or a positive whole number.
checkPatience = function(patience){

  if(!is.null(patience) && (length(patience) != 1 || !is.numeric(patience) || patience < 1 || patience != round(patience))){
    stop("patience must be NULL or a positive whole number.", call. = FALSE)
  }

  invisible(patience)
}

# Tracks path selection values as candidates are produced, so a search can stop
# once `patience` candidates in a row have not improved on the best so far.
# add() scores a new candidate and returns TRUE when the search should stop.
# With undersmooth, the best candidate is fixed when patience runs out, and the
# search continues until a candidate has twice its number of regions.
pathPatienceTracker = function(patience, type, V_Smooth, undersmooth = F){

  scoreOne = function(candidate){
    pathSelectionValues(list(candidate), type, V_Smooth)
  }

  patienceTracker(patience, scoreOne, function(values){pathBestIndex(values, type)},
                  function(candidate){length(candidate$Path)}, undersmooth)
}

# The shared patience logic for the path and trajectory searches.
# * score(candidate): the candidate's selection value (a number, or a column for BICVel).
# * best(values): the index of the best candidate so far.
# * n_regions(candidate): the candidate's number of regions.
# add() returns TRUE when the search should stop: once `patience` candidates in
# a row have not improved on the best, or, with undersmooth, once a candidate
# has at least twice the regions of the best at the moment patience ran out.
# best_index() is that best candidate (or the best overall if patience never ran out).
patienceTracker = function(patience, score, best, n_regions, undersmooth = F){

  values = NULL
  regions = c()
  fixed_best = NULL

  add = function(candidate){
    v = score(candidate)
    values <<- if(is.null(values)) v else if(is.matrix(v)) cbind(values, v) else c(values, v)
    regions <<- c(regions, n_regions(candidate))
    n = length(regions)

    if(is.null(patience)){
      return(FALSE)
    }

    if(is.null(fixed_best) && (n - best(values)) >= patience){
      fixed_best <<- best(values)
    }

    if(is.null(fixed_best)){
      return(FALSE)
    }

    !undersmooth || regions[n] >= 2*regions[fixed_best]
  }

  best_index = function(){
    if(!is.null(fixed_best)) fixed_best else if(is.null(values)) NULL else best(values)
  }

  list(add = add, values = function(){values}, best_index = best_index)
}

# The undersmoothed choice: the first candidate with at least twice the regions
# of the selected one, or the last candidate if none has that many.
undersmoothIndex = function(n_regions, best_index){

  idx = which(n_regions >= 2*n_regions[best_index])[1]

  if(is.na(idx)){
    idx = length(n_regions)
  }

  idx
}
