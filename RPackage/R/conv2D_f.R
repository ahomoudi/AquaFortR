#' @title 2D Convolution using Fortran.
#'
#' @description  This function calculates the 2D Convolution of two matrices
#'  `a` and `b` using compiled Fortran subroutine.
#'
#' @param a A matrix (2D array) of values.
#' @param b A matrix (2D array) of values.
#' @return A matrix representing the 2D cross-correlation of the input matrices.
#' @export
#' @examples
#' a <- matrix(c(1, 2, 3, 4), ncol = 2)
#' b <- matrix(c(5, 6, 7, 8), ncol = 2)
#' conv2D_f(a, b)
#' @author Ahmed Homoudi
#' @export
conv2D_f <- function(a, b) {
  stopifnot(length(dim(a)) == 2 | length(dim(b)) == 2)
  result <- .Call(
    c_conv2d_f,
    a,
    b
  )
  return(result)
}
