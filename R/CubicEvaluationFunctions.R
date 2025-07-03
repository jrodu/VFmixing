#' evaluateCubic
#'
#' @description
#' Evaluate a series of centered cubic functions at a given time.
#' \deqn{f(t) = a(t - t_0)^3 + b(t - t_0)^2 + c(t - t_0) + d}
#'
#'
#' @param t a numeric vector. The times at which the cubic functions should be evaluated.
#' @param Cubic an n x 4 numeric matrix.
#' The first columns contains the cubic coefficients (a), the second the quadratic (b), etc.
#' @param startTime a numeric scalar. The value of the centering time \eqn{t_0}.
#'
#' @returns a numeric matrix. Each row corresponds to the outputs of each cubic function's value at the time values.
#'
#' @export
#'
#' @examples
#' evaluateCubic(1,c(1,1,1,1),0)
#' evaluateCubic(0:5, matrix(c(1,1,1,1,
#'                             1,2,3,4,
#'                             2,1,6,4), nrow = 3, byrow = TRUE), 10)
evaluateCubic = function(t, Cubic, startTime){

  Cubic = matrix(Cubic, ncol = 4)

  eval = function(t){
    c(Cubic %*% c((t - startTime)^3,(t - startTime)^2,(t - startTime),1))
  }

  matrix(t(sapply(X = t, FUN = eval)), ncol = nrow(Cubic))

}

#' evaluateCubicVelocity
#'
#' @description
#' Evaluate the first derivative of a series of centered cubic functions at a given time.
#' \deqn{f(t) = a(t - t_0)^3 + b(t - t_0)^2 + c(t - t_0) + d}
#'
#' @param t a numeric vector. The times at which the cubic functions should be evaluated.
#' @param Cubic an n x 4 numeric matrix.
#' The first columns contains the cubic coefficients (a), the second the quadratic (b), etc.
#' @param startTime a numeric scalar. The value of the centering time \eqn{t_0}.
#'
#' @returns a numeric matrix. Each row corresponds to the outputs of each cubic function's first derivative at the time values.
#'
#' @export
#'
#' @examples
#' evaluateCubicVelocity(1,c(1,1,1,1),0)
#' evaluateCubicVelocity(0:5, matrix(c(1,1,1,1,
#'                                     1,2,3,4,
#'                                     2,1,6,4), nrow = 3, byrow = TRUE), 10)
evaluateCubicVelocity = function(t, Cubic, startTime){

  Cubic = matrix(Cubic, ncol = 4)

  eval = function(t){
    c(Cubic %*% c(3*(t - startTime)^2,2*(t - startTime),1,0))
  }

  matrix(t(sapply(X = t, FUN = eval)), ncol = nrow(Cubic))
}

#' evaluateCubicAcceleration
#'
#' @description
#' Evaluate the second derivative of a series of centered cubic functions at a given time.
#' \deqn{f(t) = a(t - t_0)^3 + b(t - t_0)^2 + c(t - t_0) + d}
#'
#' @param t a numeric vector. The times at which the cubic functions should be evaluated.
#' @param Cubic an n x 4 numeric matrix.
#' The first columns contains the cubic coefficients (a), the second the quadratic (b), etc.
#' @param startTime a numeric scalar. The value of the centering time \eqn{t_0}.
#'
#' @returns a numeric matrix. Each row corresponds to the outputs of each cubic function's second derivative at the time values.
#'
#' @export
#'
#' @examples
#' evaluateCubicAcceleration(1,c(1,1,1,1),0)
#' evaluateCubicAcceleration(0:5, matrix(c(1,1,1,1,
#'                                         1,2,3,4,
#'                                         2,1,6,4), nrow = 3, byrow = TRUE), 10)
evaluateCubicAcceleration = function(t, Cubic, startTime){

  Cubic = matrix(Cubic, ncol = 4)

  eval = function(t){
    c(Cubic %*% c(6*(t - startTime),2,0,0))
  }

  matrix(t(sapply(X = t, FUN = eval)), ncol = nrow(Cubic))
}

#' evaluateCubicJerk
#'
#' @description
#' Evaluate the third derivative of a series of centered cubic functions at a given time.
#' \deqn{f(t) = a(t - t_0)^3 + b(t - t_0)^2 + c(t - t_0) + d}
#'
#' @param t a numeric vector. The times at which the cubic functions should be evaluated.
#' @param Cubic an n x 4 numeric matrix.
#' The first columns contains the cubic coefficients (a), the second the quadratic (b), etc.
#' @param startTime a numeric scalar. The value of the centering time \eqn{t_0}.
#'
#' @returns a numeric matrix. Each row corresponds to the outputs of each cubic function's third derivative at the time values.
#'
#' @export
#'
#' @examples
#' evaluateCubicJerk(1,c(1,1,1,1),0)
#' evaluateCubicJerk(0:5, matrix(c(1,1,1,1,
#'                                 1,2,3,4,
#'                                 2,1,6,4), nrow = 3, byrow = TRUE), 10)
evaluateCubicJerk = function(t, Cubic, startTime){

  Cubic = matrix(Cubic, ncol = 4)

  eval = function(t){
    c(Cubic %*% c(6,0,0,0))
  }

  matrix(t(sapply(X = t, FUN = eval)), ncol = nrow(Cubic))
}

