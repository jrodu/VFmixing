#' Base Vector Fields
#'
#' @description
#' A base vector field function.
#' This is a template for creating any vector field function for using this package.
#'
#'
#' @param t a numeric scalar. The current time value.
#' @param curPos a d-dimensional numeric vector. The current position in the field.
#'
#' @returns a d x m matrix. Each column corresponds to the vector for one model.
#' 
#' @export
#' 
#' @examples
#' baseVectorFields(0, c(1,1))
baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))

  matrix(c(f1,f2), nrow = 2, byrow = F)

}

#' Flat Vector Field Function
#'
#' @description
#' Flattens the matrix output from a vector field function.
#'
#'
#' @param data_v a (d+1)-dimensional numeric vector. The time value is the first entry followed by the position.
#' @param baseVectorFields a vector field function.
#'
#' @returns a (dxm)-dimensional numeric vector.
#' The rows of the outputted matrix from baseVectorFields are concatenated.
#' For example: c(X1v1,X2v1,X1v2,X2v2)
#' 
#' @export
#'
#' @examples
#' flatVF(c(0,1,1), baseVectorFields)
#' 
flatVF = function(data_v, baseVectorFields){

  c(baseVectorFields(data_v[1],data_v[-1]))

}

#' TrajWeightedBaseVectorFields
#' 
#' @description
#' Weight a baseVectorFields function using an EstimatedTrajectory object.
#' 
#' @param t a numeric scalar. The time the vector fields should be evaluated.
#' @param curPos a numeric vector. The position in the vector field.
#' @param baseVectorFields a baseVectorFields function.
#' @param TrajList a list of numeric matrices. The list of coefficient matrices of the EstimatedTrajectory object.
#' @param TrajTimeSplits a numeric vector. The nodes of the EstimatedTrajectory object in addition with the starting and ending time.
#'
#' @returns a numeric vector. The velocity vector of the trajectory-weighted vector field at the evaluated point in time and position.
#' @export
#' 
#' @examples
#' TrajWeightedBaseVectorFields(t = 1, curPos = c(1,1), baseVectorFields, 
#'                              TrajList = ExTraj$Traj, 
#'                              TrajTimeSplits = ExTraj$TimeSplits)
#' 
TrajWeightedBaseVectorFields = function(t, curPos, baseVectorFields, TrajList, TrajTimeSplits){
  
  cur_traj_index = min(c(max(c(sum((TrajTimeSplits - t)<0),1)), length(TrajTimeSplits)-1))
  
  cur_traj = TrajList[[cur_traj_index]]
  
  cur_traj_value = evaluateCubic(t, Cubic = cur_traj, startTime = TrajTimeSplits[cur_traj_index])
  
  cur_theorVel = baseVectorFields(t, curPos)
  
  t(cur_theorVel %*% t(cur_traj_value))
  
}


