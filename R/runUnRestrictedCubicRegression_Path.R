#' runUnRestrictedCubicRegression_Path
#'
#' @param cur_Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#' @param start_time a numeric scalar. The centering time for the cubic function being estimated.
#' @param dim a numeric scalar. The dimension of the data to estimate the cubic.
#'
#' @returns a list. The elements are the estimated coefficients, the estimated value of sigma2, 
#' and the model matrix.
#' @export
#'
runUnRestrictedCubicRegression_Path = function(cur_Data, start_time, dim){
  
  cur_X = cbind((cur_Data$t - start_time)^3, 
                (cur_Data$t - start_time)^2,
                (cur_Data$t - start_time),
                1)  
  cur_Y = cur_Data[,dim+1]
  
  cur_coef = t(backsolve(qr.R(qr(cur_X)), t(qr.Q(qr(cur_X))) %*% cur_Y))
  cur_sigma2 = mean(((cur_X %*% t(cur_coef)) - cur_Y)^2)
  
  list(cur_coef,cur_sigma2,cur_X)
  
}