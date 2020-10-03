
exampledat <- list(
pat = system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid"),
mov = system.file("extdata", "Example_movement_file.csv", package="a2bcovid"),
ali = system.file("extdata", "Example_sequences.fa", package="a2bcovid"),
ward = system.file("extdata", "Example_ward_file.csv", package="a2bcovid")
)

shinyServer(function(input, output, session) {

    res <- a2bcovid(pat_file = exampledat$pat,
                    mov_file = exampledat$mov,
                    ward_file = exampledat$ward,
                    ali_file = exampledat$ali)

    output$reshead <- renderTable(
        head(res)
    )

    output$rasterPlot <- renderPlot({
        plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw")
    })

})

