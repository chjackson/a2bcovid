
pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
mov_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
ward_file <- system.file("extdata", "Example_ward_file.csv", package="a2bcovid")

library(testthat)
expect_error(
res <- a2bcovid(pat_file = pat_file, mov_file = mov_file,
             ali_file = ali_file,
             data_type = 3),   "\"ward_file\" is missing")

expect_error(
  res <- a2bcovid(pat_file = pat_file, mov_file = mov_file,
                  ali_file = ali_file, ward_file = "wibble",
                  data_type = 3),   "`wibble` not found")

res <- a2bcovid(pat_file = pat_file, mov_file = mov_file,
             ali_file = ali_file, ward_file = ward_file,
             data_type = 3)

library(stringr)
library(ggplot2)
res$consistency <- ordered(str_trim(res$consistency),
                          levels=c("Unlikely","Borderline","Consistent"))
ggplot(res, aes(y=from, x=to)) +
  geom_raster(aes(fill=consistency)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.title = element_blank())
