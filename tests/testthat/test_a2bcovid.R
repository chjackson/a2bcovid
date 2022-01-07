library(a2bcovid)

pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
hcw_loc_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
pat_loc_file <- system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid")

test_that("a2bcovid including all data sources",{
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file)
    expect_equal(res$likelihood[1], -5.89489874393478) 
    expect_equal(as.character(res$consistency[1]), "Consistent") 
    expect_equal(res$p[1], 0.248914026408368) 
})

test_that("a2bcovid with only symptom data",{
    res <- a2bcovid(pat_file = pat_file)
    expect_equal(res$likelihood[1], -2.26756960394651) 
})

test_that("a2bcovid with symptom and patient location data",{
    res <- a2bcovid(pat_file = pat_file, pat_loc_file = pat_loc_file)
    expect_equal(res$likelihood[1], -5.21925118062198) 
})

test_that("a2bcovid with symptom and staff location data",{
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file)
    expect_equal(res$likelihood[1], -2.26756960394651)
    # staff location data are uninformative in this example
})

test_that("a2bcovid with only symptom, patient location and staff location data",{
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    pat_loc_file = pat_loc_file)
    expect_equal(res$likelihood[1], -5.21925118062198) 
})

test_that("changing the virus strain",{
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file, 
                    strain = "delta")
    expect_equal(res$likelihood[1], -8.94060169337857) 
    expect_equal(res$p[1], 0.024727018479327) 
})

test_that("errors in file specification", {
    expect_error(
        res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                        ali_file = ali_file, pat_loc_file = "wibble"),
        "`wibble` not found")
    expect_error(res <- a2bcovid(pat_file = c(1,2,3)), "pat_file should be a character vector")
    expect_error(res <- a2bcovid(pat_file = pat_file, ali_file = 1, "ali_file should be"))
})

test_that("plot method runs without error",{
    expect_error({
        res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                        ali_file = ali_file, pat_loc_file = pat_loc_file)
        plot_a2bcovid(res)
        plot(res)
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=FALSE)
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=TRUE)
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw",cluster=FALSE,continuous=TRUE)
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", hi_col="purple")
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="PuRd")
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn")
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw", palette="BuGn", direction=-1)
    }, NA)
})

test_that("diagnostic output", {
    res <- a2bcovid(pat_file = pat_file, diagnostic = TRUE)
    expect_equal(res$likelihood[1], -2.26756960394651) 
})


test_that("Various alternative parameter settings", {
    res <- a2bcovid(pat_file = pat_file, chat = 0.7)
    expect_equal(res$likelihood[1], -2.26756960394651) 
    expect_equal(res$p[1], 1) 
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    evo_rate = 0.001)
    expect_equal(res$likelihood[1], -5.96038813476298)

    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    seq_noise = 0.5)
    expect_equal(res$likelihood[1], -5.98120874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    min_qual = 0.9)
    expect_equal(res$likelihood[1], -5.89489874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    max_n = 20)
    expect_equal(res$likelihood[1], -5.89489874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    pat_default = 0.5)
    expect_equal(res$likelihood[1], -5.89489874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    hcw_default = 0.5)
    expect_equal(res$likelihood[1], -5.89489874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    use_all_seqs = 1)
    expect_equal(res$likelihood[1], -5.89489874393478)
    
    res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                    ali_file = ali_file, pat_loc_file = pat_loc_file,
                    symptom_uncertainty_calc = 1)
    expect_equal(res$likelihood[1], -5.89489874393478)
})

test_that("results agree between long and wide file formats", {
    wide_file <- system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid")
    long_file <- wide_to_long(wide_file)
    res1 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                     ali_file = ali_file, pat_loc_file = long_file)
    res2 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                     ali_file = ali_file, pat_loc_file = wide_file)
    expect_equal(res1$likelihood, res2$likelihood)
    
    long_file <- system.file("extdata", "Example_pat_loc_file_long.csv", package="a2bcovid")
    res3 <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file,
                     ali_file = ali_file, pat_loc_file = long_file)
    expect_equal(res1$likelihood, res3$likelihood)
})
