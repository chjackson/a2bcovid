
pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
hcw_loc_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
pat_loc_file <- system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid")

library(testthat)

res <- a2bcovid(pat_file = pat_file)

res <- a2bcovid(pat_file = pat_file, strain="default")
res <- a2bcovid(pat_file = pat_file, strain="delta")
res <- a2bcovid(pat_file = pat_file, strain=list(pa=90, pb=0.2, po=25, smu=1.3, ssigma=0.7))
expect_error(res <- a2bcovid(pat_file = pat_file, strain=list(pa=90, pb=0.2, po=25, smu=1.3)),
             "strain\\[\\[\"ssigma\"\\]\\] should be a number")
expect_error(res <- a2bcovid(pat_file = pat_file, strain="wibble"), "`strain` should be")

res <- a2bcovid(pat_file = pat_file,
                ali_file = ali_file, pat_loc_file = pat_loc_file)

res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                ali_file = ali_file, pat_loc_file = pat_loc_file)
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw")
plot_a2bcovid(res)
plot(res)

expect_error(
  res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                  ali_file = ali_file, pat_loc_file = "wibble"),
  "`wibble` not found")
expect_error(res <- a2bcovid(pat_file = c(1,2,3)), "pat_file should be a character vector")
expect_error(res <- a2bcovid(pat_file = pat_file, ali_file = 1, "ali_file should be"))

## todo get seq or not
res$x
pdat <- approx(x = thresholds_seq$lik, y = thresholds_seq$p,
               xout=res$likelihood, yleft=1, yright=0)
do.call("cbind", pdat)

plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=FALSE)
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=TRUE)

## TODO make sure palettes match.
# Exclude diagonals first
# L > -8.15176 consistent,      L>-10.1578 borderline, else unlikely
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=FALSE)
plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=FALSE,continuous=TRUE)

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


# Different file formats

newdata_file <- "../CUH_new_patient_movement_format_edit.csv"
odfile <- newdata_to_olddata(newdata_file)
read.csv(odfile)

pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
hcw_loc_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                ali_file = ali_file, pat_loc_file = odfile)

# patients in odfile are different from patients in other files.
# Doesn't complain when different people put in different files.

wide_file <- system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid")
long_file <- wide_to_long(wide_file)
res1 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                 ali_file = ali_file, pat_loc_file = long_file)
res2 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                 ali_file = ali_file, pat_loc_file = wide_file)
expect_equal(res1, res2)

long_file <- system.file("extdata", "Example_pat_loc_file_long.csv", package="a2bcovid")
res3 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                 ali_file = ali_file, pat_loc_file = long_file)
expect_equal(res1, res3)
