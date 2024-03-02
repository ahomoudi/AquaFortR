## code to prepare `radiosonde` dataset goes here
library(readr)

radiosonde <- read_csv("data-raw/radiosonde-data.txt")[, -1] # remove row numbers
usethis::use_data(radiosonde, overwrite = TRUE)
