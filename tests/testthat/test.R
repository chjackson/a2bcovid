
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

## TODO good to bad color scheme plus colorblind safe.  inverse viridis?
## its a bit bright.


res_from <- res[!duplicated(res$from),]
hcw_from <- res_from$from_hcw[match(levels(factor(res$from)), res_from$from)]
x_cols <- ifelse(hcw_from, "red", "black")
res_to <- res[!duplicated(res$to),]
hcw_to <- res_to$to_hcw[match(levels(factor(res$to)), res_to$to)]
y_cols <- ifelse(hcw_to, "red", "black")

ggplot(res, aes(y=from, x=to)) +
  geom_raster(aes(fill=consistency)) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust=0.5, colour = x_cols),
        axis.text.y = ggtext::element_markdown(colour = y_cols),
        legend.title = element_blank())


