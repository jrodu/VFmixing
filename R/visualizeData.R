#' visualizeData
#'
#' @param Data a numeric data frame. The first column is the time and the remaining columns are the positions.
#'
#' @returns a ggplot2 object.
#' @export
#'
visualizeData = function(Data){
  ggplot2::ggplot(data = Data, ggplot2::aes(x = X1, y = X2)) + ggplot2::geom_point()
}