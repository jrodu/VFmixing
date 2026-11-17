#' constructUnboundedCubicSplineTransitionMatrix_Traj
#'
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#' @param V_Smooth a boolean. TRUE if the trajectory should be estimated using C2 smoothness. FALSE if C1.
#'
#' @returns a numeric matrix. Transforms the free trajectory coefficients into the full coefficients.
#' @export
#'
constructUnboundedCubicSplineTransitionMatrix_Traj = function(time_splits, n_models, V_Smooth = T){

  Path_TransMat = constructUnboundedCubicSplineTransitionMatrix_Path(time_splits, V_Smooth)

  # Each region's rows are a block-diagonal copy of the path transition rows,
  # one block per model: filled into a preallocated matrix (binding row by row
  # copies the whole matrix for every region)
  n_regions = length(time_splits) - 1
  n_free = ncol(Path_TransMat)
  Traj_TransMat = matrix(0, 4*n_models*n_regions, n_models*n_free)

  for(i in seq_len(n_regions)){
    block = Path_TransMat[1:4 + (i-1)*4, , drop = FALSE]
    for(m in seq_len(n_models)){
      Traj_TransMat[(i-1)*4*n_models + (m-1)*4 + 1:4, (m-1)*n_free + seq_len(n_free)] = block
    }
  }

  Traj_TransMat

}
