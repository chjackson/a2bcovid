
library(whoinfected)
library(tidyverse)

setwd("/Users/chris/Dropbox/work/covid/cill/whoinfected/tests/local")

res <- mainR(pat_file = "cluster_D_genetic_temporal_data_NoPII_20200807.csv",
      mov_file = "hcw_movement_20200811_D_NoPII.csv",
      data_type = 3)

res$consistency <- ordered(str_trim(res$consistency),
                          levels=c("Unlikely","Borderline","Consistent"))
ggplot(res, aes(y=from, x=to)) +
  geom_raster(aes(fill=consistency)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.title = element_blank())

OurApp()
