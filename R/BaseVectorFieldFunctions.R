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

  f1 = c(curPos[2],-1*curPos[1])
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
#' baseVectorFields = function(t, curPos){
#'
#' f1 = c(curPos[2],-1*curPos[1])
#' f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
#'
#' matrix(c(f1,f2), nrow = 2, byrow = FALSE)
#'
#'}
#' 
#' flatVF(c(0,1,1), baseVectorFields)
flatVF = function(data_v, baseVectorFields){

  c(baseVectorFields(data_v[1],data_v[-1]))

}




