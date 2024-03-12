conv2D_f0 <- function(a, b) {
  dyn.load("fortran/conv2D.so")

  # the full convolution matrix
  conv_row <- nrow(a) + nrow(b) - 1
  conv_col <- ncol(a) + ncol(b) - 1
  conv <- matrix(1:c(conv_row * conv_col),
    byrow = FALSE, ncol = conv_col
  )

  result <- .Fortran("conv2d_f",
    m = as.integer(dim(a)[1]),
    n = as.integer(dim(a)[2]),
    p = as.integer(dim(b)[1]),
    q = as.integer(dim(b)[2]),
    k = as.integer(conv_row),
    l = as.integer(conv_col),
    a = as.double(a),
    b = as.double(b),
    conv = as.double(conv)
  )$conv

  return(matrix(result, nrow = conv_row, ncol = conv_col))
}
