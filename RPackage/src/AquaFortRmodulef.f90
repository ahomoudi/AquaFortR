module AquaFortRmodule
   use, intrinsic :: iso_c_binding
   implicit none

   public :: xcorr2d_f
   public :: conv2d_f
contains

   ! Cross-Correlation
   subroutine xcorr2d_f(m, n, p, q, k, l, a, b, cc) bind(C, name = "xcorr2d_f_")
      implicit none
      integer(kind=c_int), intent(in), value                :: m, n, p, q, k, l
      real(kind=c_double), intent(in), dimension(m, n)      :: a
      real(kind=c_double), intent(in), dimension(p, q)      :: b
      real(kind=c_double), intent(out), dimension(k, l)     :: cc
      !     dummy vars
      integer(kind=c_int)                                   :: min_row_shift, min_col_shift
      integer(kind=c_int)                                   :: max_row_shift, max_col_shift
      integer(kind=c_int)                                   :: rows_padded, cols_padded
      integer(kind=c_int)                                   :: icol, irow, icc, icol2, irow2
      real(kind=c_double), allocatable, dimension(:, :)     :: padded_a, padded_b

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
            icc = irow + ((icol - 1)*k)
            icol2 = icol + q - 1
            irow2 = irow + p - 1
            padded_b(irow:irow2, icol:icol2) = b
            cc(irow, icol) = sum(padded_a*padded_b)
            padded_b = 0.0
         end do
      end do
   end subroutine xcorr2d_f

   ! Convolution 

   subroutine conv2d_f(m, n, p, q, k, l, a, b, conv) bind(C, name = "conv2d_f_")
      implicit none
      integer(kind=c_int), intent(in), value                :: m, n, p, q, k, l
      real(kind=c_double), intent(in), dimension(m, n)      :: a
      real(kind=c_double), intent(in), dimension(p, q)      :: b
      real(kind=c_double), intent(out), dimension(k, l)     :: conv
      !     dummy vars
      integer(kind=c_int)                                   :: min_row_shift, min_col_shift, i, j
      integer(kind=c_int)                                   :: max_row_shift, max_col_shift
      integer(kind=c_int)                                   :: rows_padded, cols_padded
      integer(kind=c_int)                                   :: icol, irow, iconv, icol2, irow2
      real(kind=c_double), allocatable, dimension(:, :)     :: padded_a, padded_b
      real(kind=c_double), dimension(p, q)                  :: b_flipped

      !     obtain possible shfits
      min_row_shift = -1*(p - 1)
      max_row_shift = m - 1
      min_col_shift = -1*(q - 1)
      max_col_shift = n - 1

      !   Flip the kernel i.e. B
      b_flipped = 0.0
      do i = 1, p
         do j = 1, q
            b_flipped(p - i + 1, q - j + 1) = b(i, j)
         end do
      end do

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
            padded_b(irow:irow2, icol:icol2) = b_flipped
            conv(irow, icol) = sum(padded_a*padded_b)
            padded_b = 0.0
         end do
      end do
   end subroutine conv2d_f

   ! CAPE

end module AquaFortRmodule
