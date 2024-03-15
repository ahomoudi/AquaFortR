conv2D_f0 <- function(a, b) {
  require(dotCall64)
  dyn.load("AquaFortR_Codes/conv2D.so")
  
  m <- nrow(a)
  n <- ncol(b)
  
  p <- nrow(b)
  q <- ncol(b)
  # the full convolution matrix
  conv_row <- m + p - 1
  conv_col <- n + q - 1
  conv <- matrix(0,
                 ncol = conv_col,
                 nrow = conv_row)
  
  conv <- .C64("conv2d_f",
               SIGNATURE = c(rep("integer",6),
                             rep("double",3)),
               INTENT = c(rep("r",8), "rw"),
               m, n, p, q,
               k = conv_row,
               l = conv_col,
               a = a, b = b,
               conv = conv)$conv
  
  return(conv)
}