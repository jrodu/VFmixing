#' Generate C1 Smoothness
#' 
#' @description
#' Generate C1 Smoothness for a cubic function given a previous cubic to connect.
#' 
#'
#' @param prevCubic an nx4 numeric matrix.
#'  Each row corresponds the coefficients of the cubic functions for which the current cubic will be connected. 
#' @param curCubic an nx2 numeric matrix. 
#' Each row corresponds the cubic and quadratic coefficients of the cubic functions currently being connected.
#' @param prevStartTime a numeric scalar. The centering time of the previous cubic functions. 
#' @param curStartTime a numeric scalar. The centering time of the current cubic functions.
#' @param meetTime a numeric scalar. The time at which the previous and current cubic functions meet.
#'
#' @returns an nx4 numeric matrix. 
#' Each row corresponds the coefficients of the current cubic functions with C1 smoothness to the previous cubic functions.
#' @export
#'
#' @examples
#' 
#' prevCubic <- matrix(c(1,2,3,4), nrow = 1, byrow = TRUE)
#' curCubic <- matrix(c(2,1/2), nrow = 1, byrow = TRUE)
#' prevStartTime <- 0
#' curStartTime <- 1
#' meetTime <- 1
#' 
#' generateC1Smoothness(prevCubic, curCubic, prevStartTime, curStartTime, meetTime)
#' 
generateC1Smoothness = function(prevCubic, curCubic, prevStartTime, curStartTime, meetTime){
  
  c = prevCubic %*% c(3*(meetTime - prevStartTime)^2, 2*(meetTime - prevStartTime), 1, 0) - 
    curCubic[,c(1,2)] %*% c(3*(meetTime - curStartTime)^2, 2*(meetTime - curStartTime))
  
  d = prevCubic %*% c((meetTime - prevStartTime)^3, (meetTime - prevStartTime)^2, (meetTime - prevStartTime), 1) - 
    curCubic[,c(1,2)] %*% c((meetTime - curStartTime)^3, (meetTime - curStartTime)^2) - c*(meetTime - curStartTime)
  
  
  curCubic = cbind(matrix(curCubic[,c(1,2)], ncol = 2), c, d)
  colnames(curCubic) = c("Beta1","Beta2","Beta3","Beta4")
  curCubic
}
