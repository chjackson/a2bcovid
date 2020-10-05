options(shiny.maxRequestSize=100*1024^2)

txtfiles_accepted <-   c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                         'text/plain', '.csv', '.tsv')

exampledat <- list(
    pat = system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid"),
    mov = system.file("extdata", "Example_movement_file.csv", package="a2bcovid"),
    ali = system.file("extdata", "Example_sequences.fa", package="a2bcovid"),
    ward = system.file("extdata", "Example_ward_file.csv", package="a2bcovid")
)

exampleres <- a2bcovid(pat_file = exampledat$pat, mov_file = exampledat$mov,
                       ward_file = exampledat$ward, ali_file = exampledat$ali)

vals <- reactiveValues(data = "default")

shinyServer(function(input, output, session) {

    observe({
        input$pat; input$mov; input$ali; input$ward
        vals$data <- "user"
    })

    get_userres <- reactive({
        dat <- exampledat
        if (!is.null(input$pat)) dat$pat <- input$pat$datapath
        if (!is.null(input$mov)) dat$mov <- input$mov$datapath
        if (!is.null(input$ali)) dat$ali <- input$ali$datapath
        if (!is.null(input$ward)) dat$ward <- input$ward$datapath
        a2bcovid(pat_file = dat$pat, mov_file = dat$mov,
                 ward_file = dat$ward, ali_file = dat$ali)
    })

    get_res <- reactive({
        get_userres()
        switch(vals$data,
               default = exampleres,
               user = get_userres())
    })

    output$reshead <- renderTable({
        head(get_res())
    })

    output$rasterPlot <- renderPlot({
        plot_a2bcovid(get_res(), hi_from="from_hcw", hi_to="to_hcw")
    })

    output$patInput <- renderUI({
        input$reset
        vals$data <- "default"
        fileInput('pat',
                  'Patient file in CSV format. [ Info about format here ].',
                  accept = txtfiles_accepted)
    })
    output$movInput <- renderUI({
        input$reset
        vals$data <- "default"
        fileInput('mov',
                  'Movement file in CSV format. [ Info about format here ].',
                  accept = txtfiles_accepted)
    })
    output$wardInput <- renderUI({
        input$reset
        vals$data <- "default"
        fileInput('ward',
                  'Ward file in CSV format. [ Info about format here ].',
                  accept = txtfiles_accepted)
    })
    output$aliInput <- renderUI({
        input$reset
        vals$data <- "default"
        fileInput('ali',
                  'Sequence file in FASTA format. [ Info about format here ].',
                  accept = ".fa")
    })
})

