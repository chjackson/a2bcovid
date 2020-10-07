
shinyUI(fluidPage(
  tags$head(
    tags$style(".resetButton {float:left; font-size:12px;}")
  ),

  titlePanel("Infection Pathways from Sequence Data"),

  sidebarLayout(

      sidebarPanel(
        h3("Upload data"),
        uiOutput('patInput'),
        uiOutput('aliInput'),
        checkboxInput("use_ali",label="Use genomic data",value=TRUE),
        uiOutput('wardInput'),
        checkboxInput("use_ward",label="Use patient location data",value=TRUE),
        uiOutput('movInput'),
        checkboxInput("use_mov",label="Use staff location data",value=FALSE),
        HTML("<button id='reset' class='action-button resetButton'>Reset to example data</button><p>")
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
