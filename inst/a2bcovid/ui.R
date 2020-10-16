data_warning <- "Please ensure the data contain no information that might disclose the identities of people, and that you are permitted to upload the data, unless you are running this app locally from R"

shinyUI(fluidPage(
  tags$head(
    tags$style(".resetButton {float:left; font-size:12px;}")
  ),

  titlePanel("a2bcovid"),

  p(code("a2bcovid"), "is a tool to determine the likelihood that an infection was transmitted between pairs of individuals, and identify clusters of infections"),
  p("Upload your data on the left, and specify which of the uplodaded datasets will be used in the calculation. Results will appear on the right."),

  sidebarLayout(

      sidebarPanel(
        bsTooltip("patInput", data_warning, placement = "bottom", trigger = "hover", options=NULL),
        bsTooltip("aliInput", data_warning, placement = "bottom", trigger = "hover", options=NULL),
        bsTooltip("wardInput", data_warning, placement = "bottom", trigger = "hover", options=NULL),
        bsTooltip("movnput", data_warning, placement = "bottom", trigger = "hover", options=NULL),
        h3("Upload data"),
        uiOutput('patInput'),
        checkboxInput("use_ali",label="Use genomic data",value=TRUE),
        uiOutput('aliInput'),
        checkboxInput("use_ward",label="Use patient location data",value=TRUE),
        uiOutput('wardInput'),
        checkboxInput("use_mov",label="Use staff location data",value=FALSE),
        uiOutput('movInput'),
        HTML("<button id='reset' class='action-button resetButton'>Reset to example data</button><p>")
      ),

      mainPanel(
        tabsetPanel(type="tabs",
                    tabPanel("Output",
                             p('The graph indicates the likelihood that infection was transmitted between each pair of individuals'),
                             textOutput('whichdata'),
                             plotOutput(outputId = "rasterPlot"),
                             checkboxInput("cluster",label="Cluster",value=TRUE),
                             checkboxInput("continuous",label="Colour gradient",value=FALSE),
                             bsCollapse(id = "collapseExample", open = "Panel 2",
                                        bsCollapsePanel("Advanced settings", p("Not implemented yet"),
                                                        p("Here you will be able to change parameters used in the algorithm"),
                                                        numericInput("seq_noise","Parameter governing sequencing noise",
                                                                     0.772469, min = NA, max = NA,  step = NA, width = NULL),
                                                        numericInput("evo_rate","Virus evolution rate (nucleotide substitutions per locus per year)",
                                                                     0.0008,  min = NA, max = NA,  step = NA, width = NULL),
                                                        numericInput("min_qual","Fraction of genome covered by the sequence",
                                                                     0.8,  min = NA, max = NA,  step = NA, width = NULL),
                                                        numericInput("max_n","Maximum number of ambiguous nucleotides tolerated",
                                                                     10, min = NA, max = NA,  step = NA, width = NULL),
                                                        style = "default")
                             )
                    ),

                    tabPanel("Help on upload format",
                             h4("Times of symptom onset"),
                             p("A .csv file with columns in this order. Column names don't matter."),
                             tags$ol(tags$li("Individual ID"),
                                     tags$li("Symptom onset date, as dd/mm/YYYY"),
                                     tags$li("Onset date source. Should take the value 1, 2 or 3:",
                                             tags$ol(tags$li("Date of onset known"),
                                                     tags$li("Asymptomatic infection"),
                                                     tags$li("Date of onset missing or unknown"))),
                                     tags$li("Infection type. Should take the value 1, 2 or 3:",
                                             tags$ol(
                                               tags$li("Patient, infected in the community TODO def"),
                                               tags$li("Patient, not infected in the community"),
                                               tags$li("Healthcare worker"))),
                                     tags$li("Sequence ID, matching ID in sequence file"),
                                     tags$li("Date of sample collection, as dd/mm/YYYY"),
                                     tags$li("Columns 7 AND ONWARDS? are ignored")),
                             p("For example:"),
                             uiOutput('pat_data_table'),

                             h4("Genome sequence file"),
                             p("A FASTA format file containing all required sequences,
                               identified by the sequence ID in the symptom onset data file."),

                             h4("Patient location data"),
                             p("A .csv file with columns in this order."),
                             tags$ol(tags$li("Individual ID (matching those in symptom onset file)"),
                                     tags$li("Ward ID"),
                                     tags$li("Infection type. This must have one of two values:",
                                             tags$ul(tags$li("patient: patient"),
                                                     tags$li("hcw: healthcare worker"))),
                                     tags$li("Availability of data e.g. \"patient_moves_available\". ????"),
                                     tags$li("Columns 4 and onwards give data on patient location",
                                             tags$br(),
                                             "These are given in groups of three columns:",
                                             tags$ol(tags$li("Name of the location"),
                                                     tags$li("Start date for individual being in that location"),
                                                     tags$li("End date for individual being in that location"))),
                                     tags$li("In practice only the first column, and columns from 5 onwards are used. [????? we shouldn't ask for things that are not needed]")
                             ),
                             uiOutput('ward_data_table'),

                             h4("Staff location data"),
                             p("A .csv file with columns in this order"),
                             tags$ol(tags$li("Individual ID (matching those in symptom onset file)"),
                                     tags$li("Ward ID"),
                                     tags$li("Columns 3 and onwards represent dates.",
                                             tags$br(),
                                             "The name of the column should be a date in the format dd.mm.YYYY [ NOTE inconsistent with other date format ] ",
                                             tags$br(),
                                             "The values in the column should be either Y or N:",
                                             tags$ul(
                                               tags$li("Y: Healthcare worker was on the ward on this date"),
                                               tags$li("N: Healthcare worker was not on the ward on this date")
                                             )
                                     )
                             ),
                             uiOutput('mov_data_table')
                    ),
                    tabPanel("Example analysis",includeHTML("example.html"))
        )
      )
  )

))
