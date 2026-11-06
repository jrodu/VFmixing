#' visualizeData
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#' @param x_var a string. The variable name to plot on the x-axis. Columns are named t, X1, X2, ...
#' @param y_var a string. The variable name to plot on the y-axis.
#'
#' @returns a ggplot2 object.
#' @export
#' 
#' @examples
#' visualizeData(ExData)
#' visualizeData(ExData, x_var = "t", y_var = "X1")
#' 
visualizeData = function(Data, x_var = "X1", y_var = "X2"){
  
  Data = data.frame(Data)
  names(Data) = c("t", stringr::str_c("X", 1:(ncol(Data)-1)))
  
  if(!all(c(x_var, y_var) %in% names(Data))){
    stop("x_var and y_var must be among: ", paste(names(Data), collapse = ", "), call. = FALSE)
  }
  
  ggplot2::ggplot(data = Data, ggplot2::aes(x = !!dplyr::sym(x_var), y = !!dplyr::sym(y_var))) + ggplot2::geom_point()
}
