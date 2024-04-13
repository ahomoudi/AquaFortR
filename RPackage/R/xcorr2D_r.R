#' @title 2D cross-correlation using R.
#'
#' @description  This function calculates the 2D cross-correlation of two matrices
#'  `a` and `b` using compiled Fortran subroutine.
#'
#' @param a A matrix (2D array) of values.
#' @param b A matrix (2D array) of values.
#' @return A matrix representing the 2D cross-correlation of the input matrices.
#' @export
#' @examples
#' a <- matrix(c(1, 2, 3, 4), ncol = 2)
#' b <- matrix(c(5, 6, 7, 8), ncol = 2)
#' xcorr2D_r(a, b)
#' @author Ahmed Homoudi
#' @export
xcorr2D_r <- function(a, b) {
  stopifnot(length(dim(a)) == 2 | length(dim(b)) == 2)
  # the full CC matrix
  cc_row <- nrow(a) + nrow(b) - 1
  cc_col <- ncol(a) + ncol(b) - 1
  cc <- matrix(1:c(cc_row * cc_col), byrow = FALSE, ncol = cc_col)

  # obtain possible shifts
  min_row_shift <- -(nrow(b) - 1)
  max_row_shift <- (nrow(a) - 1)
  min_col_shift <- -(ncol(b) - 1)
  max_col_shift <- (ncol(a) - 1)

  # Padded matrix
  rows_padded <- abs(min_row_shift) + nrow(a) + abs(max_row_shift)
  cols_padded <- abs(min_col_shift) + ncol(a) + abs(max_col_shift)
  # a
  padded_a <- matrix(0, nrow = rows_padded, ncol = cols_padded)
  padded_a[
    (abs(min_row_shift) + 1):(abs(min_row_shift) + nrow(a)),
    (abs(min_col_shift) + 1):(abs(min_col_shift) + ncol(a))
  ] <- a

  for (icol in 1:cc_col) {
    for (irow in 1:cc_row) {
      icc <- irow + ((icol - 1) * cc_row)
      cols <- (icol):(icol + ncol(b) - 1)
      rows <- (irow):(irow + nrow(b) - 1)
      # b
      padded_b <- array(0, dim = c(rows_padded, cols_padded))
      padded_b[rows, cols] <- b

      cc[irow, icol] <- sum(padded_a * padded_b)
    }
  }

  return(cc)
}
