library(data.table)
library(tidyverse)

data<-fread("/media/SSD2/DATA/DWD/sekundenwerte_aero_04928_20210101_20211231_hist/produkt_sec_aero_20210101_20211231_04928.txt",
            nThread = 6)

data <- data %>%
  mutate(BEZUGSDATUM_SYNOP = strptime(BEZUGSDATUM_SYNOP, format ="%Y%m%d%H"))


example_data <- data %>%
  filter(BEZUGSDATUM_SYNOP=="2021060112")

fwrite(example_data, "/media/SSD2/DATA/DWD/sekundenwerte_aero_04928_20210101_20211231_hist/radiosonde-data.txt")
ggplot(example_data)+
  geom_point(aes(x = AE_GB_POS,
                 y = MESSZEITPUNKT,
                 color = AE_RF))

library(car)
scatter3d(x = example_data$AE_GL_POS, y = example_data$AE_GB_POS, z = example_data$MESSZEITPUNKT,
          surface=FALSE)


as.POSIXct(x = as.vector(),
           format ="%Y%m%d%H",
           tz = "CET")

#data$BEZUGSDATUM_SYNOP<- as.Date(data$BEZUGSDATUM_SYNOP, format = "%Y%m%d%H")

library(stringr)
data$BEZUGSDATUM_SYNOP1<- as.POSIXct(data$BEZUGSDATUM_SYNOP,
                                  "%Y%m%d%H",
                                  tz = "CET")
hist(data$BEZUGSDATUM_SYNOP1, "months")
library(ggplot2)


ggplot(data_June)+
  geom_point(aes(x = BEZUGSDATUM_SYNOP1,
                       y = MESSZEITPUNKT,
                 color = AE_TT))

data_June<- data[which(month(data$BEZUGSDATUM_SYNOP1)==6),]

final_data<-data_June%>%
  dplyr::filter(BEZUGSDATUM_SYNOP1 %in% sample(unique(data_June$BEZUGSDATUM_SYNOP1), 5))

unique(final_data$BEZUGSDATUM_SYNOP)
fwrite(final_data, "/media/SSD2/DATA/DWD/sekundenwerte_aero_04928_20210101_20211231_hist/produkt_sec_aero_20210601_20210630_04928.txt")
