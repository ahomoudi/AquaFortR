module AquaFortRmodule
   use, intrinsic :: iso_c_binding
   implicit none

   public :: xcorr2d_f
   public :: conv2d_f
   public :: cape_f
contains

   ! Cross-Correlation
   subroutine xcorr2d_f(m, n, p, q, k, l, a, b, cc) bind(C, name="xcorr2d_f_")
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

   subroutine conv2d_f(m, n, p, q, k, l, a, b, conv) bind(C, name="conv2d_f_")
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
   subroutine cape_f(t_parcel, dwpt_parcel, mr_parcel, &
                     nlevel, p_profile, t_profile, mr_profile, &
                     vtc, nresult, result) bind(C, name="cape_f_")
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
      implicit none
      real(c_double), intent(in)                           :: t_parcel, dwpt_parcel, mr_parcel
      integer(c_int), value                                :: nlevel, nresult
      real(c_double), intent(in), dimension(nlevel)        :: p_profile, t_profile, mr_profile
      integer(c_int), value                                :: vtc
      ! result CAPE, CIN, p_LCL, p_LFC
      real(c_double), intent(out), dimension(nresult)  :: result

      ! Constants
      real(c_double), parameter:: Rd = 287.058 ! gas constant for dry air [J * kg**-1 * K**-1]
      real(c_double), parameter:: Rv = 461.5   ! gas constant for water vapour [J * kg**-1 * K**-1]

      ! Dummy variables
      integer(c_int)             :: ilevel
      real(c_double)             :: dp, dt, pF, rF
      real(c_double)             :: t_parcel_buoyancy, t_env_buoyancy, t_parcel_tmp
      real(c_double)             :: CIN_add, CAPE_add, rF1, p_EL
      logical(c_bool)            :: ice = .false.
      !real(c_double), parameter  :: NaN = TRANSFER((/Z'00000000', Z'7FF80000'/), 1.0_8)
      real(c_double)             :: NaN

      ! Starting values
      NaN = ieee_value(1., ieee_quiet_nan)
      t_parcel_tmp = t_parcel
      ilevel = 1
      result(1) = 0.0 ! CAPE
      result(2) = 0.0 ! CIN
      result(3) = NaN ! LCL
      result(4) = NaN ! LFC
      CIN_add = 0.0
      CAPE_add = 0.0

      ! Find LCL & Estimate CIN
      do
         if ((t_parcel_tmp .le. dwpt_parcel) .or. (ilevel .ge. nlevel)) exit
         !change in pressure
         dp = p_profile(ilevel + 1) - p_profile(ilevel)
         ! change in temperature
         dt = calc_dtdp_dry(t_parcel_tmp, p_profile(ilevel))*dp
         ! new parcel temperature
         t_parcel_tmp = t_parcel_tmp + dt
         ! liquid mixing ratio
         rF = mr_parcel
         if (vtc .eq. 1) then
            ! Früh and Wirth Eq. 12: virtual temperature coorection to calc CAPE
            t_parcel_buoyancy = t_parcel_tmp*(1.0 + (Rv/Rd)*rF)/(1.0 + rF)
            t_env_buoyancy = t_profile(ilevel) + (1.0 + (Rv/Rd)*mr_profile(ilevel))/ &
                             (1.0 + mr_profile(ilevel))
         else
            t_parcel_buoyancy = t_parcel_tmp
            t_env_buoyancy = t_profile(ilevel)
         end if
         ! Accumlate CIN
         CIN_add = -1.0*Rd*((t_parcel_buoyancy - t_env_buoyancy)* &
                            (log(p_profile(ilevel)) - log(p_profile(ilevel + 1))))
         result(2) = result(2) + CIN_add
         ilevel = ilevel + 1
         !print*, ilevel, t_parcel_tmp
      end do
      result(3) = p_profile(ilevel) ! LCL
      ! find LFC
      if (t_parcel_tmp .le. t_profile(ilevel)) then
         do
            if ((t_parcel_tmp .ge. t_profile(ilevel)) .or. (ilevel .ge. nlevel)) exit
            dp = p_profile(ilevel + 1) - p_profile(ilevel)
            !
            pF = saturation_vapour_pressure(t_parcel_tmp, ice)
            ! Saturation mixing ratio
            rF = calc_rF1(pF, p_profile(ilevel) - pF)
            !
            dt = calc_dtdp_wet(t_parcel_tmp, p_profile(ilevel), rF)*dp
            ! new parcel temperature
            t_parcel_tmp = t_parcel_tmp + dt
            if (vtc .eq. 1) then
               ! Früh and Wirth Eq. 12: virtual temperature coorection to calc CAPE
               t_parcel_buoyancy = t_parcel_tmp*(1.0 + (Rv/Rd)*rF)/(1.0 + rF)
               t_env_buoyancy = t_profile(ilevel) + (1.0 + (Rv/Rd)*mr_profile(ilevel))/ &
                                (1.0 + mr_profile(ilevel))
            else
               t_parcel_buoyancy = t_parcel_tmp
               t_env_buoyancy = t_profile(ilevel)
            end if
            ! Accumlate CIN
            CIN_add = -1.0*Rd*((t_parcel_buoyancy - t_env_buoyancy)* &
                               (log(p_profile(ilevel)) - log(p_profile(ilevel + 1))))
            result(2) = result(2) + CIN_add
            ilevel = ilevel + 1
            !print*, ilevel, t_parcel_tmp
         end do
         if (ilevel .lt. (nlevel - 1)) then
            result(4) = p_profile(ilevel)
         else
            result(4) = NaN
         end if

      else
         result(4) = p_profile(ilevel)
      end if

      ! find EL and estimate CAPE
      do
         if ((t_parcel_tmp .le. t_profile(ilevel)) .or. (ilevel .ge. nlevel)) exit
         !
         pF = saturation_vapour_pressure(t_parcel_tmp, ice)
         ! Saturation mixing ratio
         rF = calc_rF1(pF, p_profile(ilevel) - pF)
         if (vtc .eq. 1) then
            ! Früh and Wirth Eq. 12: virtual temperature coorection to calc CAPE
            t_parcel_buoyancy = t_parcel_tmp*(1.0 + (Rv/Rd)*rF)/(1.0 + rF)
            t_env_buoyancy = t_profile(ilevel) + (1.0 + (Rv/Rd)*mr_profile(ilevel))/ &
                             (1.0 + mr_profile(ilevel))
         else
            t_parcel_buoyancy = t_parcel_tmp
            t_env_buoyancy = t_profile(ilevel)
         end if
         CAPE_add = -1.0*Rd*((t_parcel_buoyancy - t_env_buoyancy) &
                             *(log(p_profile(ilevel)) - log(p_profile(ilevel + 1))))

         result(1) = result(1) + CAPE_add
         dp = p_profile(ilevel + 1) - p_profile(ilevel)
         rF1 = 0.0
         dt = calc_dtdp_wet(t_parcel_tmp, p_profile(ilevel), rF1)*dp
         t_parcel_tmp = t_parcel_tmp + dt
         ilevel = ilevel + 1
         !print*, ilevel, t_parcel_tmp
      end do

      result(1) = abs(result(1))
      !
      if (ilevel .eq. nlevel) then
         p_EL = NaN
      else
         p_EL = p_profile(ilevel)
      end if
      ! Calculate upper part above EL
      do
         if (ilevel .ge. nlevel) exit
         rF1 = 0.0
         dp = p_profile(ilevel + 1) - p_profile(ilevel)
         dt = calc_dtdp_wet(t_parcel_tmp, p_profile(ilevel), rF1)*dp
         t_parcel_tmp = t_parcel_tmp + dt
         ilevel = ilevel + 1
         !print*, ilevel, t_parcel_tmp
      end do
      !print*, p_EL
   end subroutine cape_f

   ! utilities
   function specific_heat_dry_air(temperature) result(C_pd)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      real(c_double), intent(in) :: temperature
      real(c_double)             :: t_temp, C_pd
      ! temperature should be: -40°C < T < 40^C
      ! output is in [J kg^-1 C^-1]
      t_temp = temperature - 273.15
      C_pd = 1005.60 + 0.017211*t_temp + 0.000392*t_temp**2.0

   end function specific_heat_dry_air

   function specific_heat_water_vapour(temperature) result(C_pv)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      real(c_double), intent(in) :: temperature
      real(c_double)             :: t_temp, C_pv

      ! Reid, R.C., J.M. Prausnitz, and B.E. Poling (1987)
      ! The Properties of Gases and Liquids.  4th ed.  McGraw-Hill, 741 pp.
      ! output is in J kg^-1 K^-1

      t_temp = temperature - 273.15
      c_pv = 1858.0 + 3.820*10.0**(-1.0)*t_temp + 4.220*10.0**(-4.0)*t_temp**2.0 - &
             1.996*10.0**(-7.0)*temperature**3.0

   end function

   function specific_heat_liquid_water(temperature) result(C_pw)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      real(c_double), intent(in) :: temperature
      real(c_double)             :: t_temp, C_pw
      ! temperature [K]
      ! output is in J kg^-1 K^-1

      t_temp = temperature - 273.15
      c_pw = 4217.4 - 3.720283*t_temp + 0.1412855*t_temp**2.0 - 2.654387*10.0**(-3.0)*t_temp**3.0 &
             + 2.093236*10.0**(-5.0)*t_temp**(4.0)

   end function

   function saturation_vapour_pressure(temperature, ice_in) result(es)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      ! calculates the saturation vapour pressure in hPa using the Clausius-Claperon equation
      real(c_double), intent(in):: temperature ! [K]
      logical(c_bool), intent(in), optional:: ice_in

      ! keyword ice, indicates if even in case of temperatures lower than 273.15 K es is calculated with
      ! respect to liquid water (then ice must not been set)
      ! output is in hPa
      ! written by K.Barfus 12/2009

      logical:: ice

      real(c_double), parameter:: e0 = 0.611 ! [kPa]
      real(c_double), parameter:: T0 = 273.15 ! [K]
      real(c_double), parameter:: Rv = 461.0 ! [J K**-1 kg**-1] gas constant for water vapour
      real(c_double):: L, es

      if (present(ice_in)) then
         ice = ice_in
      else
         ice = .false.
      end if

      if (ice .eqv. (.true.)) then
         if (temperature .gt. 273.15) then  ! water
            L = 2.5*10.0**6.0 ! J kg**-1
         else
            L = 2.83*10.0**6.0  ! J kg**-1
         end if
      else
         L = 2.5*10.0**6.0 ! J kg**-1
      end if

      if (temperature .gt. 0) then
         es = e0*exp((L/Rv)*(1.0/T0 - 1.0/temperature))
         es = es*10.0
      else
         es = 0.0
      end if

   end function saturation_vapour_pressure

   function calc_rF1(pF1, p0) result(res) ! Frueh and Wirth, Eq. 4
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none

      real(c_double), intent(in):: pF1  ! saturation vapour pressure [hPa]
      real(c_double), intent(in):: p0   ! partial pressure of dry air [hPa]

      real(c_double), parameter:: Rd = 287.058 ! gas constant for dry air [J * kg**-1 * K**-1]
      real(c_double), parameter:: Rv = 461.5   ! gas constant for water vapour [J * kg**-1 * K**-1]

      real(c_double)           :: res

      res = (Rd*pF1)/(Rv*p0)

   end function

   function latent_heat_gas_to_liquid(temperature) result(res)
      ! latent heat of condensation due to Rogers and Yau in J/kg
      ! valid for 248.15 K < T < 313.15 K
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      real(c_double), intent(in):: temperature ! temperature [K]

      real(c_double):: t_temp, latent_heat, res

      t_temp = temperature - 273.15
      latent_heat = 2500.8 - 2.36*t_temp + 0.0016*t_temp**2.0 - 0.00006*t_temp**3.0
      res = latent_heat*1000.0

      ! alternative approach
      ! calculates the latent heat of condensation (gas -> liquid) due to
      ! Fleagle, R.G. and J.A. Businger, (1980)
      ! An Introduction to Atmospheric Physics.  2d ed.  Academic Press, 432 pp.
      ! input
      ! T in K
      ! output in J kg^-1 K^-1
      ! t_temp = T - 273.15
      ! Lv = (25.00 - 0.02274 * t_temp) * 10.0^5.0

   end function

   function calc_dtdp_dry(temperature, pressure) result(dtdp)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      real(c_double), intent(in):: temperature ! temperature [K]
      real(c_double), intent(in):: pressure ! pressure [hPa]
      ! output is:
      ! [K/hPa]

      real(c_double), parameter:: Rd = 287.058 ! gas constant for dry air [J * kg**-1 * K**-1]
      real(c_double):: cp0, dtdp

      cp0 = specific_heat_dry_air(temperature)
      dtdp = (temperature*Rd)/(pressure*cp0)

   end function

   function calc_dtdp_wet(temperature, pressure, rF) result(dtdp)
      ! Written by Klemens Barfus. Modified by Ahmed Homoudi
      implicit none
      ! not applying mixed-phase model !
      real(c_double), intent(in):: temperature  ! temperature [K]
      real(c_double), intent(in):: pressure  ! pressure [hPa]
      real(c_double), intent(in), optional:: rF ! liquid mixing ratio [kg/kg]
      ! rF is liquid mixing ratio <- here 0.0 because of an irreversible process
      ! output is:
      ! [K/hPa]
      real(c_double), parameter::  Rd = 287.058 ! gas constant for dry air [J * kg**-1 * K**-1]
      real(c_double), parameter::  Rv = 461.5   ! gas constant for water vapour [J * kg**-1 * K**-1]
      real(c_double):: pF1, rF1, lF1, LLF1
      real(c_double):: cp0, cp1, cp2, Cp, v, dtdp, p0

      pF1 = saturation_vapour_pressure(temperature) ! hPa
      p0 = pressure - pF1                      ! hPa
      rF1 = calc_rF1(pF1, p0)  ! saturation mixing ratio in g/g
      lF1 = latent_heat_gas_to_liquid(temperature) ! J/kg
      LLF1 = pF1*lF1                   ! hPa * (J/kg)
      cp0 = specific_heat_dry_air(temperature)     ! J/(kg*K)
      cp1 = specific_heat_water_vapour(temperature)  ! J/(kg*K)
      cp2 = specific_heat_liquid_water(temperature)  ! J/(kg*K)
      Cp = cp0 + cp1*rF1 + cp2*rF      ! J/(kg*K)
      v = (rF1*lF1)/pF1*(1.0 + (Rv/Rd)*rF1)*(LLF1/(Rv*temperature**2.0))
      dtdp = ((rF1*Rv*temperature/pF1)*(1.0 + (rF1*lF1/(Rd*temperature))))/(Cp + v)

   end function
end module AquaFortRmodule
