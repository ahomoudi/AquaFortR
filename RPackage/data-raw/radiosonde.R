## code to prepare `radiosonde` dataset goes here
library(readr)

radiosonde <- read_csv("data-raw/radiosonde-data.txt",
  col_types = cols(BEZUGSDATUM_SYNOP = col_datetime(format = "%Y%m%d%H"))
)
usethis::use_data(radiosonde, overwrite = TRUE)
