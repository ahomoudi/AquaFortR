xcorr2D_f0 <- function(a, b) {
  # Please adjust the path to your setup.
  dyn.load("AquaFortR_Codes/xcorr2D.so")

  # the full CC matrix
  cc_row <- nrow(a) + nrow(b) - 1
  cc_col <- ncol(a) + ncol(b) - 1
  cc <- matrix(1:c(cc_row * cc_col), byrow = FALSE, ncol = cc_col)

  result <- .Fortran("xcorr2d_f",
    m = as.integer(dim(a)[1]),
    n = as.integer(dim(a)[2]),
    p = as.integer(dim(b)[1]),
    q = as.integer(dim(b)[2]),
    k = as.integer(cc_row),
    l = as.integer(cc_row),
    a = as.double(a),
    b = as.double(b),
    cc = as.double(cc)
  )$cc

  return(matrix(result, nrow = nrow(cc), ncol = ncol(cc)))
}
