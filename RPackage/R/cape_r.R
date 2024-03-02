#' @title CAPE using R
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
#' cape_r(t_profile[1], dwpt_profile[1], mr_profile[1],
#'   p_profile, t_profile, mr_profile,
#'   vtc = TRUE
#' )
#' @author Klemens Barfus (Original in Fortran), Ahmed Homoudi (Translation to R)
#' @useDynLib AquaFortR
#' @export
cape_r <- function(t_parcel, dwpt_parcel, mr_parcel,
                   p_profile, t_profile, mr_profile,
                   vtc = TRUE) {
  # Constants
  Rd <- 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  Rv <- 461.5 # gas constant for water vapour [J * kg**-1 * K**-1]


  t_parcel_tmp <- t_parcel
  ilevel <- 1
  nlevel <- length(p_profile)
  CAPE <- 0.0
  CIN <- 0.0
  p_LCL <- NaN
  p_LFC <- NaN

  # LCL and CIN
  while ((t_parcel_tmp > dwpt_parcel) & ilevel < nlevel) {
    # change in pressure
    dp <- p_profile[ilevel + 1] - p_profile[ilevel]
    #  change in temperature
    dt <- calc_dtdp_dry(t_parcel_tmp, p_profile[ilevel]) * dp
    # new parcel temperature
    t_parcel_tmp <- t_parcel_tmp + dt
    # liquid mixing ratio
    rF <- mr_parcel
    if (vtc) {
      # Früh and Wirth Eq. 12: virtual temperature correction to calc CAPE
      t_parcel_buoyancy <- t_parcel_tmp * (1.0 + (Rv / Rd) * rF) / (1.0 + rF)
      t_env_buoyancy <- t_profile[ilevel] + (1.0 + (Rv / Rd) * mr_profile[ilevel]) /
        (1.0 + mr_profile[ilevel])
    } else {
      t_parcel_buoyancy <- t_parcel_tmp
      t_env_buoyancy <- t_profile[ilevel]
    }
    CIN_add <- -1.0 * Rd * ((t_parcel_buoyancy - t_env_buoyancy) *
      (log(p_profile[ilevel]) - log(p_profile[ilevel + 1])))

    CIN <- CIN + CIN_add
    ilevel <- ilevel + 1
  }
  p_LCL <- p_profile[ilevel]

  # LFC
  if (t_parcel_tmp < t_profile[ilevel]) {
    while (t_parcel_tmp < t_profile[ilevel] & ilevel < nlevel) {
      # change in pressure
      dp <- p_profile[ilevel + 1] - p_profile[ilevel]
      #
      pF <- saturation_vapour_pressure(t_parcel_tmp)
      rF <- calc_rF1(pF, p_profile[ilevel] - pF)
      #  change in temperature
      dt <- calc_dtdp_wet(t_parcel_tmp, p_profile[ilevel], rF) * dp
      # new parcel temperature
      t_parcel_tmp <- t_parcel_tmp + dt
      if (vtc) {
        # Früh and Wirth Eq. 12: virtual temperature correction to calc CAPE
        t_parcel_buoyancy <- t_parcel_tmp * (1.0 + (Rv / Rd) * rF) / (1.0 + rF)
        t_env_buoyancy <- t_profile[ilevel] + (1.0 + (Rv / Rd) * mr_profile[ilevel]) /
          (1.0 + mr_profile[ilevel])
      } else {
        t_parcel_buoyancy <- t_parcel_tmp
        t_env_buoyancy <- t_profile[ilevel]
      }
      CIN_add <- -1.0 * Rd * ((t_parcel_buoyancy - t_env_buoyancy) *
        (log(p_profile[ilevel]) - log(p_profile[ilevel + 1])))

      CIN <- CIN + CIN_add
      ilevel <- ilevel + 1
    }
    if (ilevel < (nlevel - 1)) {
      p_LFC <- p_profile[ilevel]
    } else {
      p_LFC <- NaN
    }
  } else {
    p_LFC <- p_profile[ilevel]
  }

  # EL and CAPE
  while (t_parcel_tmp > t_profile[ilevel] & ilevel < nlevel) {
    pF <- saturation_vapour_pressure(t_parcel_tmp)
    rF <- calc_rF1(pF, p_profile[ilevel] - pF)
    if (vtc) {
      # Früh and Wirth Eq. 12: virtual temperature correction to calc CAPE
      t_parcel_buoyancy <- t_parcel_tmp * (1.0 + (Rv / Rd) * rF) / (1.0 + rF)
      t_env_buoyancy <- t_profile[ilevel] + (1.0 + (Rv / Rd) * mr_profile[ilevel]) /
        (1.0 + mr_profile[ilevel])
    } else {
      t_parcel_buoyancy <- t_parcel_tmp
      t_env_buoyancy <- t_profile[ilevel]
    }
    CAPE_add <- -1.0 * Rd * ((t_parcel_buoyancy - t_env_buoyancy) *
      (log(p_profile[ilevel]) - log(p_profile[ilevel + 1])))
    CAPE <- CAPE + CAPE_add
    dp <- p_profile[ilevel + 1] - p_profile[ilevel]
    rF1 <- 0
    dt <- calc_dtdp_wet(t_parcel_tmp, p_profile[ilevel], rF1) * dp
    t_parcel_tmp <- t_parcel_tmp + dt
    ilevel <- ilevel + 1
  }

  CAPE <- abs(CAPE)
  if (ilevel == nlevel) {
    p_EL <- NaN
  } else {
    p_EL <- p_profile[ilevel]
  }
  while (ilevel < nlevel) {
    rF1 <- 0
    dp <- p_profile[ilevel + 1] - p_profile[ilevel]
    dt <- calc_dtdp_wet(t_parcel_tmp, p_profile[ilevel], rF1) * dp
    t_parcel_tmp <- t_parcel_tmp + dt
    ilevel <- ilevel + 1
  }

  # END
  result <- c(CAPE, CIN, p_LCL, p_LFC)
  names(result) <- c("CAPE", "CIN", "p_LCL", "p_LFC")
  return(result)
}

