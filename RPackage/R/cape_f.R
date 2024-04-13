#' @title CAPE using Fortran
#'
#' @description  The function calculates the Convective Available Potential
#'  Energy (CAPE) using R. It returns a vector of CAPE, convective inhibition
#'  (CIN), pressure at lifted condensation level (LCL),and level of free
#'  convection (LFC).
#'
#' @param t_parcel Parcel temperature [K].
#' @param dwpt_parcel Parcel dew point temperature [K].
#' @param mr_parcel Parcel mixing ratio [kg/kg].
#' @param p_profile Pressure profile from sounding or modelling [hPa].
#' @param t_profile Temperature profile from sounding or modelling [K].
#' @param mr_profile Mixing ratio profile from sounding or modelling [kg/kg].
#' @param vtc logical refers is virtual temperature correction due to Doswell and
#' Rasmussen (1994).
#' @return A vector containing CAPE, CIN, p_LCL, p_LFC
#' @examples
#' data("radiosonde")
#' t_profile <- radiosonde$temp + 273.15 # K
#' dwpt_profile <- radiosonde$dpt + 273.15 # K
#' p_profile <- radiosonde$pressure # hPa
#' # Mixing ratio
#' mr_profile <- mixing_ratio_from_dewpoint(dwpt_profile, p_profile)
#' cape_f(t_profile[1], dwpt_profile[1], mr_profile[1],
#'   p_profile, t_profile, mr_profile,
#'   vtc = TRUE
#' )
#' @author Klemens Barfus (Original in Fortran), Ahmed Homoudi (Integration in R)
#' @export
cape_f <- function(t_parcel, dwpt_parcel, mr_parcel,
                   p_profile, t_profile, mr_profile,
                   vtc = TRUE) {

  vtc <-as.integer(vtc)

  result <- .Call(
    c_cape_f,
    t_parcel,
    dwpt_parcel,
    mr_parcel,
    p_profile,
    t_profile,
    mr_profile,
    vtc
  )

  names(result) <- c("CAPE", "CIN", "p_LCL", "p_LFC")
  return(result)
}


