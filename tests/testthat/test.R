
pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
mov_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
ward_file <- system.file("extdata", "Example_ward_file.csv", package="a2bcovid")

library(testthat)

expect_error(
  res <- a2bcovid(pat_file = pat_file, mov_file = mov_file,
                  ali_file = ali_file, ward_file = "wibble"),   "`wibble` not found")

res <- a2bcovid(pat_file = pat_file, mov_file = mov_file,
             ali_file = ali_file, ward_file = ward_file)

plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", cluster=FALSE)
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", hi_col="purple")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="PuRd")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn")
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn", direction=-1)

# Available palettes.  Select direction = 1 or -1 to reverse the colours
# Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn,
# PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd

set.seed(1)
res2 <- res[sample(1:nrow(res), replace=FALSE),]
res2$from <- factor(as.character(res2$from), levels=unique(as.character(res2$from)))
res2$to <- factor(as.character(res2$to), levels=unique(as.character(res2$to)))

plot_a2bcovid(res2, hi_from="from_hcw", hi_to="to_hcw")
plot_a2bcovid(res2, hi_from="from_hcw", hi_to="to_hcw", cluster=FALSE)

if (0){

  ## Run this to test the web app locally
  a2bcovid_app()  # from installed package

  shiny::runApp("inst/a2bcovid")  # directly from the R code

}
