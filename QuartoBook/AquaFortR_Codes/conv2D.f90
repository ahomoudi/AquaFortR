subroutine conv2d_f(m, n, p, q, k, l, a, b, conv)
   implicit none
   integer                                :: m, n, p, q, k, l, i, j
   double precision, dimension(m, n)      :: a
   double precision, dimension(p, q)      :: b
   double precision, dimension(k, l)      :: conv
   !     dummy vars
   integer                               :: min_row_shift, min_col_shift
   integer                               :: max_row_shift, max_col_shift
   integer                               :: rows_padded, cols_padded
   integer                               :: icol, irow, iconv, icol2, irow2
   real, allocatable, dimension(:, :)    :: padded_a, padded_b

   !     obtain possible shfits
   min_row_shift = -1*(p - 1)
   max_row_shift = m - 1
   min_col_shift = -1*(q - 1)
   max_col_shift = n - 1


   !   Padded arrray
   rows_padded = abs(min_row_shift) + m + abs(max_row_shift)
   cols_padded = abs(min_col_shift) + n + abs(max_col_shift)
   !    A
   allocate (padded_a(rows_padded, cols_padded))
   padded_a = 0.0
   padded_a((abs(min_row_shift) + 1):(abs(min_row_shift) + m), &
            (abs(min_col_shift) + 1):(abs(min_col_shift) + n)) = a

   !    B
   allocate (padded_b(rows_padded, cols_padded))
   padded_b = 0.0
   do icol = 1, l
      do irow = 1, k
         iconv = irow + ((icol - 1)*k)
         icol2 = icol + q - 1
         irow2 = irow + p - 1
         padded_b(irow:irow2, icol:icol2) = b(p:1:-1,q:1:-1)
         conv(irow, icol) = sum(padded_a*padded_b)
         padded_b = 0.0
      end do
   end do
end subroutine conv2d_f
