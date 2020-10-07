options(shiny.maxRequestSize=100*1024^2)

txtfiles_accepted <-   c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                         'text/plain', '.csv', '.tsv')

exampledat <- list(
    pat = system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid"),
    mov = system.file("extdata", "Example_movement_file.csv", package="a2bcovid"),
    ali = system.file("extdata", "Example_sequences.fa", package="a2bcovid"),
    ward = system.file("extdata", "Example_ward_file.csv", package="a2bcovid")
)

vals <- reactiveValues(data = "default",
                       using_gen  = TRUE,
                       using_ward = TRUE,
                       using_mov = FALSE)

shinyServer(function(input, output, session) {

    observe({
        input$pat; input$mov; input$ali; input$ward
        vals$using_gen  <- input$dataType %in% c("times_gen","times_gen_pat","times_gen_pat_staff")
        vals$using_ward <- input$dataType %in% c("times_gen_pat","times_gen_pat_staff")
        vals$using_mov  <- input$dataType %in% c("times_gen_pat_staff")
        if (!is.null(input$pat))  vals$data <- "user"
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
        if (!is.null(input$mov)) dat$mov <- input$mov$datapath
        if (!is.null(input$ali)) dat$ali <- input$ali$datapath
        if (!is.null(input$ward)) dat$ward <- input$ward$datapath
        a2bcovid(pat_file = dat$pat, mov_file = dat$mov,
                 ward_file = dat$ward, ali_file = dat$ali,
                 data_type = match(input$dataType, data_choices) - 1)
    })

    get_exampleres <- reactive({
        a2bcovid(pat_file = exampledat$pat, mov_file = exampledat$mov,
                 ward_file = exampledat$ward, ali_file = exampledat$ali,
                 data_type =  match(input$dataType, data_choices) - 1)
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
        input$dataType
        fileInput('ali',
                  'Genome sequence file (FASTA)',
                  accept = ".fa")
    })
    output$wardInput <- renderUI({
        input$reset
        input$dataType
        fileInput('ward',
                  'Patient location data (CSV)',
                  accept = txtfiles_accepted)
    })
    output$movInput <- renderUI({
        input$reset
        input$dataType
        fileInput('mov',
                  'Staff location data (CSV)',
                  accept = txtfiles_accepted)
    })

    output$whichdata <- renderText({
        patf <-  "Patient"
        alif <-  if (vals$using_gen) " | Sequence" else ""
        wardf <-  if (vals$using_ward) " | Patient location" else ""
        movf <-  if (vals$using_mov) " | Staff location" else ""
        if (vals$data=="default") {
            outstr <- sprintf("Using built-in example data:\n%s%s%s%s", patf, alif, wardf, movf)
        } else if (vals$data == "user"){
            outstr <- sprintf("Using user data:\n%s%s%s%s", patf, alif, wardf, movf)
        }
        else outstr <- "ERROR 1"
        outstr
    })
})

