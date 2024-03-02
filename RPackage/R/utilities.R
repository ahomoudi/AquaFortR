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
  return(es*10)
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
#' @title Mixing ratio
#' @description
#' Estimation of mixing ratio from pressure [hPa] and dewpoint temperature [K].
#' @param dewpoint temperature [K]
#' @param pressure  in [hPa]
#' @author Ahmed Homoudi
#' @export
mixing_ratio_from_dewpoint <- function(dewpoint, pressure) {
  eipslon <- 0.622 # g g**-1 or kg kg **-1
  e_in <- saturation_vapour_pressure(dewpoint)
  r <- (eipslon * e_in) / (pressure - e_in)
  return(r)
}
# dimensionless
mixing_ratio <- function(P_in, e_in) {
  eipslon <- 0.622 # g g**-1 or kg kg **-1
  r <- (eipslon * e_in) / (P_in - e_in)
  return(r)
}
