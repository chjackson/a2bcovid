options(shiny.maxRequestSize=100*1024^2)

txtfiles_accepted <-   c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                         'text/plain', '.csv', '.tsv')

exampledat <- list(
    pat = system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid"),
    hcwloc = system.file("extdata", "Example_hcw_loc_file.csv", package="a2bcovid"),
    ali = system.file("extdata", "Example_sequences.fa", package="a2bcovid"),
    patloc = system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid"),
    patloclong = system.file("extdata", "Example_pat_loc_file_long.csv", package="a2bcovid")
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
        if (!is.null(input$ali) && input$use_ali)
            dat$ali <- input$ali$datapath
        else {
            dat$ali <- ""
            updateCheckboxInput(session, "use_ali", value=FALSE)
        }
        if (!is.null(input$patloc) && input$use_patloc)
            dat$patloc <- input$patloc$datapath
        else {
            updateCheckboxInput(session, "use_patloc", value=FALSE)
            dat$patloc <- ""
        }
        if (!is.null(input$hcwloc) && input$use_hcwloc)
            dat$hcwloc <- input$hcwloc$datapath
        else {
            dat$hcwloc <- ""
            updateCheckboxInput(session, "use_hcwloc", value=FALSE)
        }
        a2bcovid(pat_file = dat$pat, hcw_loc_file = dat$hcwloc,
                 pat_loc_file = dat$patloc, ali_file = dat$ali,
                 seq_noise = input$seq_noise, 
                 max_n = input$max_n, min_qual = input$min_qual)
    })

    get_exampleres <- reactive({
        a2bcovid(pat_file = exampledat$pat,
                 hcw_loc_file = if (input$use_hcwloc) exampledat$hcwloc else "",
                 pat_loc_file = if (input$use_patloc) exampledat$patloc else  "",
                 ali_file = if (input$use_ali) exampledat$ali else "",
                 seq_noise = input$seq_noise, 
                 max_n = input$max_n, min_qual = input$min_qual,
                 strain = input$virus_strain)
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
        plot_a2bcovid(get_res(), cluster=input$cluster,
                      hi_from="from_hcw", hi_to="to_hcw",
                      continuous=(input$colour_scale=="continuous")
                      )
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
    output$patLocInput <- renderUI({
        input$reset
        fileInput('patloc',
                  'Patient location data (CSV)',
                  accept = txtfiles_accepted)
    })
    output$hcwLocInput <- renderUI({
        input$reset
        fileInput('hcwloc',
                  'Staff location data (CSV)',
                  accept = txtfiles_accepted)
    })

    output$whichdata <- renderText({
        patf <-  "Patient"
        alif <-  if (input$use_ali) " | Sequence" else ""
        patlocf <-  if (input$use_patloc) " | Patient location" else ""
        hcwlocf <-  if (input$use_hcwloc) " | Staff location" else ""
        if (vals$data=="default") {
            outstr <- sprintf("Using built-in example data:\n%s%s%s%s", patf, alif, patlocf, hcwlocf)
        } else if (vals$data == "user"){
            outstr <- sprintf("Using user data:\n%s%s%s%s", patf, alif, patlocf, hcwlocf)
        }
        else outstr <- "ERROR 1"
        outstr
    })

    output$pat_data_table <- renderTable({
        read.csv(exampledat$pat)[1:2,]
    })
    output$patloc_data_table <- renderTable({
        read.csv(exampledat$patloc)[1:2,]
    })
    output$patloclong_data_table <- renderTable({
        read.csv(exampledat$patloclong)[1:2,]
    })
    output$hcwloc_data_table <- renderTable({
        read.csv(exampledat$hcwloc, check.names = FALSE)[1:2,]
    })

    observe({
        input$reset_pars
        for (i in c("seq_noise","max_n","min_qual"))
            updateNumericInput(session, i, value=defaultpars[[i]])
    })
})

