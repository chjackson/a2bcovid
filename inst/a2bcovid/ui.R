
shinyUI(fluidPage(
  tags$head(
    tags$style(".resetButton {float:left; font-size:12px;}")
  ),

  titlePanel("Infection Pathways from Sequence Data"),

  sidebarLayout(

      sidebarPanel(
        uiOutput('patInput'),
        uiOutput('movInput'),
        uiOutput('aliInput'),
        uiOutput('wardInput'),
        HTML("<button id='reset' class='action-button resetButton'>Reset to example data</button>"),
        selectInput('dataType',
                    label = "Data to include",
                    choices = data_choices,
                    selected = data_choices[3])
      ),

      mainPanel(
          h4("Output here"),
          tableOutput(outputId = "reshead"),
          plotOutput(outputId = "rasterPlot")
      )
  )

))
