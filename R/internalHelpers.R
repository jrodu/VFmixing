# Internal helpers shared by the path and trajectory functions.

# Smoothness of an estimated path spline, inferred from the number of free
# parameters per dimension: 4 + (R - 1) for C2, 4 + 2(R - 1) for C1, where R is
# the number of regions. ByFours paths store no free parameters and are always C1.
inferPathVSmooth = function(EstPath, default = FALSE){

  n_regions = length(EstPath$Path)

  if(is.null(EstPath$Free_Parameters) || n_regions < 2){
    return(default)
  }

  n_free = ncol(EstPath$Free_Parameters)

  if(n_free == 4 + (n_regions - 1)){
    TRUE
  } else if(n_free == 4 + 2*(n_regions - 1)){
    FALSE
  } else{
    default
  }

}

# Number of free spline parameters per dimension (or model).
splineFreeParams = function(n_regions, V_Smooth){

  4 + V_Smooth*(n_regions - 1) + (!V_Smooth)*2*(n_regions - 1)

}

# A single error summary for an EstimatedTrajectory. ByFours trajectories store a
# matrix of squared velocity residuals (time in the first column); OneShot
# trajectories store a vector of integrated squared errors per region.
trajErrorStat = function(cur_traj){

  Error = cur_traj$Error

  if(is.null(dim(Error))){
    sum(Error)/(max(cur_traj$TimeSplits) - min(cur_traj$TimeSplits))
  } else{
    stats::median(rowSums(as.matrix(Error[,-1, drop = FALSE])))
  }

}

# Report the selected model in the same way for every selection function.
reportBestRegions = function(n_regions){

  print(stringr::str_c("Best number of regions was ", n_regions, "."))

}
