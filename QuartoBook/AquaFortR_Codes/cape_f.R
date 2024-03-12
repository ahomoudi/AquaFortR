cape_f0 <- function(t_parcel, dwpt_parcel, mr_parcel,
                    p_profile, t_profile, mr_profile,
                    vtc = TRUE) {
  dyn.load("fortran/cape_f.so")

  nlevel <- length(p_profile)
  nresult <- 4

  result <- .Fortran("cape_f",
    t_parcel = as.double(t_parcel),
    dwpt_parcel = as.double(dwpt_parcel),
    mr_parcel = as.double(mr_parcel),
    nlevel = as.integer(nlevel),
    p_profile = as.double(p_profile),
    t_profile = as.double(t_profile),
    mr_profile = as.double(mr_profile),
    vtc = as.integer(vtc),
    nresult = as.integer(nresult),
    result = as.double(rep(0, nresult))
  )$result

  # END
  names(result) <- c("CAPE", "CIN", "p_LCL", "p_LFC")
  return(result)
}
