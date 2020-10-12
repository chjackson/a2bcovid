
pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
mov_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
ward_file <- system.file("extdata", "Example_ward_file.csv", package="a2bcovid")

pat_file

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

plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", hi_col="purple")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="PuRd")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn", direction=-1)

# Available palettes.  Select direction = 1 or -1 to reverse the colours
# Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn,
# PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd

if (0){

  ## Run this to test the web app locally
  a2bcovid_app()

}
