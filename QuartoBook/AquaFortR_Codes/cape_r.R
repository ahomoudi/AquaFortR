cape_r0 <- function(t_parcel, dwpt_parcel, mr_parcel,
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

specific_heat_dry_air <- function(temperature) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # temperature [K]
  t_temp <- temperature - 273.15
  C_pd <- 1005.60 + (0.017211 * t_temp) + (0.000392 * t_temp**2.0)
  return(C_pd)
}

specific_heat_water_vapour <- function(temperature) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # temperature [K]
  t_temp <- temperature - 273.15
  c_pv <- 1858.0 + (3.820 * 10.0**(-1.0) * t_temp) +
    (4.220 * 10.0**(-4.0) * t_temp**2.0) -
    (1.996 * 10.0**(-7.0) * temperature**3.0)
  return(c_pv)
}

specific_heat_liquid_water <- function(temperature) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # temperature [K]
  t_temp <- temperature - 273.15
  c_pw <- 4217.4 - (3.720283 * t_temp) +
    (0.1412855 * t_temp**2.0) -
    (2.654387 * 10.0**(-3.0) * t_temp**3.0) +
    (2.093236 * 10.0**(-5.0) * t_temp**(4.0))
  return(c_pw)
}
#' @title Saturation vapour pressure
#' @description
#' Estimation of Saturation vapour pressure [hPa] from temperature [k].
#' @param temperature in [k]
#' @param ice TRUE or FALSE, if to consider ice state or not.
#' @author Klemens Barfus
#' @export
saturation_vapour_pressure <- function(temperature, ice = FALSE) {
  e0 <- 0.611 # kPa
  Rv <- 461.0 # J K**-1 kg**-1
  T0 <- 273.15 # K
  Lv <- 2.50e6 # J kg **-1
  Ld <- 2.83e6 # J kg **-1

  if (ice) {
    if (temperature > T0) {
      L <- Lv
    } else {
      L <- Ld
    }
  } else {
    L <- Lv
  }

  es <- ifelse(temperature > 0,
    e0 * exp((L / Rv) * ((1 / T0) - (1 / temperature))),
    0
  )
  return(es * 10)
}

calc_rF1 <- function(pF1, p0) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # Constants
  Rd <- 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  Rv <- 461.5 # gas constant for water vapour [J * kg**-1 * K**-1]
  return((Rd * pF1) / (Rv * p0))
}
latent_heat_gas_to_liquid <- function(temperature) {
  # temperature [K]
  t_temp <- temperature - 273.15
  latent_heat <- 2500.8 - 2.36 * t_temp + 0.0016 * t_temp**2.0 - 0.00006 * t_temp**3.0
  return(latent_heat * 1000.0)
}
calc_dtdp_dry <- function(temperature, pressure) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # temperature [K]
  # pressure [hPa]

  Rd <- 287.058 # gas constant for dry air [J * kg**-1 * K**-1]

  cp0 <- specific_heat_dry_air(temperature)
  dtdp <- (temperature * Rd) / (pressure * cp0)
  return(dtdp)
}
calc_dtdp_wet <- function(temperature, pressure, rF) {
  # Written by Klemens Barfus. Translated to R by Ahmed Homoudi
  # Constants
  Rd <- 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  Rv <- 461.5 # gas constant for water vapour [J * kg**-1 * K**-1]

  pF1 <- saturation_vapour_pressure(temperature)
  p0 <- pressure - pF1
  rF1 <- calc_rF1(pF1, p0)
  lF1 <- latent_heat_gas_to_liquid(temperature)
  LLF1 <- pF1 * lF1
  cp0 <- specific_heat_dry_air(temperature)
  cp1 <- specific_heat_water_vapour(temperature)
  cp2 <- specific_heat_liquid_water(temperature)
  Cp <- cp0 + cp1 * rF1 + cp2 * rF
  v <- (rF1 * lF1) / pF1 * (1.0 + (Rv / Rd) * rF1) * (LLF1 / (Rv * temperature**2.0))
  dtdp <- ((rF1 * Rv * temperature / pF1) * (1.0 + (rF1 * lF1 / (Rd * temperature)))) / (Cp + v)
  return(dtdp)
}
