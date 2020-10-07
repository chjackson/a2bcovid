
shinyUI(fluidPage(
  tags$head(
    tags$style(".resetButton {float:left; font-size:12px;}")
  ),

  titlePanel("Infection Pathways from Sequence Data"),

  sidebarLayout(

      sidebarPanel(
        h3("Upload data"),
        uiOutput('patInput'),
        conditionalPanel(
          condition = "(input.dataType == 'times_gen') || (input.dataType == 'times_gen_pat') || (input.dataType == 'times_gen_pat_staff')",
          uiOutput('aliInput')
        ),
        conditionalPanel(
          condition = "(input.dataType == 'times_gen_pat') || (input.dataType == 'times_gen_pat_staff')",
          uiOutput('wardInput')
        ),
        conditionalPanel(
          condition = "input.dataType == 'times_gen_pat_staff'",
          uiOutput('movInput')
        ),
        HTML("<button id='reset' class='action-button resetButton'>Reset to example data</button><p>"),
        selectInput('dataType',
                    label = "Data to include",
                    choices = data_choices,
                    selected = data_choices[3])
      ),

      mainPanel(
        tabsetPanel(type="tabs",
                    tabPanel("Output",
                             h4("Output here"),
                             tableOutput(outputId = "reshead"),
                             textOutput('whichdata'),
                             plotOutput(outputId = "rasterPlot"),
                             checkboxInput("cluster",label="Cluster",value=TRUE)
                    ),
                    tabPanel("Help on upload format",
                             p("TODO"),
                             h4("Times of symptom onset"),
                             h4("Genome sequences"),
                             h4("Patient location data"),
                             h4("Staff location data")
                    )
        )
      )
  )

))
