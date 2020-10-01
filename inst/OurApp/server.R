
exampledir <- "/Users/chris/work/covid/cill/a2bcovid/tests/local"
exampledat <- list(pat="cluster_D_genetic_temporal_data_NoPII_20200807.csv",
                   mov="hcw_movement_20200811_D_NoPII.csv",
                   ward="ward_movement_network_edit_anonymised_20200811_NoPII.csv",
                   ali="Seqs_editN20_manali_plus.fa")
exampledat <- lapply(exampledat, function(x)file.path(exampledir, x))

shinyServer(function(input, output, session) {

    res <- mainR(pat_file = exampledat$pat,
                 mov_file = exampledat$mov,
                 ward_file = exampledat$ward,
                 ali_file = exampledat$ali)
    res$consistency <- ordered(str_trim(res$consistency),
                               levels=c("Unlikely","Borderline","Consistent"))

    output$reshead <- renderTable(
        head(res)
    )

    output$rasterPlot <- renderPlot({
        ggplot(res, aes(y=from, x=to)) +
            geom_raster(aes(fill=consistency)) +
            theme(axis.text.x = element_text(angle = 90, vjust=0.5),
                  legend.title = element_blank())
    })

})

