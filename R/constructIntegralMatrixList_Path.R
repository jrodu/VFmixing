#' constructIntegralMatrixList_Path
#'
#' @param cur_integral_list a list of matrices. Each corresponds with a block of the previous overall integral matrix.
#' @param time_splits a numeric vector. The nodes of the current spline with the starting and ending times. 
#' @param worst_split a numeric scalar. The region to be split.
#'
#' @returns a list of matrices. 
#' @export
#'
constructIntegralMatrixList_Path = function(cur_integral_list, time_splits, worst_split){
  
  integrateMatrixElement = function(lower, upper, index_sum){
    
    (upper-lower)^(9-index_sum)/(9-index_sum)
    
  }
  
  first_elements = unlist(apply(matrix(2:8), MARGIN = 1, FUN = integrateMatrixElement, lower = time_splits[worst_split], upper = time_splits[worst_split+1]))
  second_elements = unlist(apply(matrix(2:8), MARGIN = 1, FUN = integrateMatrixElement, lower = time_splits[worst_split+1], upper = time_splits[worst_split+2]))
  
  first_elements = first_elements[c(1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,7)]
  second_elements = second_elements[c(1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,7)]
  
  first_block = matrix(first_elements, nrow = 4, byrow = T)
  second_block = matrix(second_elements, nrow = 4, byrow = T)
  
  cur_integral_list[[worst_split]] = first_block
  cur_integral_list = append(cur_integral_list, list(second_block), after = worst_split)
  
  cur_integral_list
  
}
