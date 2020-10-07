options(shiny.maxRequestSize=100*1024^2)

txtfiles_accepted <-   c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                         'text/plain', '.csv', '.tsv')

exampledat <- list(
    pat = system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid"),
    mov = system.file("extdata", "Example_movement_file.csv", package="a2bcovid"),
    ali = system.file("extdata", "Example_sequences.fa", package="a2bcovid"),
    ward = system.file("extdata", "Example_ward_file.csv", package="a2bcovid")
)

vals <- reactiveValues(data = "default")

shinyServer(function(input, output, session) {

    observe({
        if (!is.null(input$pat))
            vals$data <- "user"
    })

    observe({
        input$reset
        vals$data <- "default"
    })

    get_userres <- reactive({
        dat <- exampledat
        if (!is.null(input$pat)) {
            dat$pat <- input$pat$datapath
            vals$data <- "user"
            print("Using user data")
        }
        ## TODO Switch off input$use_mov if no data supplied.
        if (!is.null(input$ali) && input$use_ali)
            dat$ali <- input$ali$datapath
        else {
            dat$ali <- ""
            updateCheckboxInput(session, "use_ali", value=FALSE)
        }
        if (!is.null(input$ward) && input$use_ward)
            dat$ward <- input$ward$datapath
        else {
            updateCheckboxInput(session, "use_ward", value=FALSE)
            dat$ward <- ""
        }
        if (!is.null(input$mov) && input$use_mov)
            dat$mov <- input$mov$datapath
        else {
            dat$mov <- ""
            updateCheckboxInput(session, "use_mov", value=FALSE)
        }
        a2bcovid(pat_file = dat$pat, mov_file = dat$mov,
                 ward_file = dat$ward, ali_file = dat$ali)
    })

    get_exampleres <- reactive({
        a2bcovid(pat_file = exampledat$pat,
                 mov_file = if (input$use_mov) exampledat$mov else "",
                 ward_file = if (input$use_ward) exampledat$ward else  "",
                 ali_file = if (input$use_ali) exampledat$ali else "")
    })

    get_res <- reactive({
        switch(vals$data,
               default = get_exampleres(),
               user = get_userres())
    })

    output$reshead <- renderTable({
        head(get_res())
    })

    output$rasterPlot <- renderPlot({
        plot_a2bcovid(get_res(), cluster=input$cluster, hi_from="from_hcw", hi_to="to_hcw")
    })

    output$patInput <- renderUI({
        input$reset
        fileInput('pat',
                  'Times of symptom onset (CSV)',
                  accept = txtfiles_accepted)
    })
    output$aliInput <- renderUI({
        input$reset
        fileInput('ali',
                  'Genome sequence file (FASTA)',
                  accept = ".fa")
    })
    output$wardInput <- renderUI({
        input$reset
        fileInput('ward',
                  'Patient location data (CSV)',
                  accept = txtfiles_accepted)
    })
    output$movInput <- renderUI({
        input$reset
        fileInput('mov',
                  'Staff location data (CSV)',
                  accept = txtfiles_accepted)
    })

    output$whichdata <- renderText({
        patf <-  "Patient"
        alif <-  if (input$use_ali) " | Sequence" else ""
        wardf <-  if (input$use_ward) " | Patient location" else ""
        movf <-  if (input$use_mov) " | Staff location" else ""
        if (vals$data=="default") {
            outstr <- sprintf("Using built-in example data:\n%s%s%s%s", patf, alif, wardf, movf)
        } else if (vals$data == "user"){
            outstr <- sprintf("Using user data:\n%s%s%s%s", patf, alif, wardf, movf)
        }
        else outstr <- "ERROR 1"
        outstr
    })
})

