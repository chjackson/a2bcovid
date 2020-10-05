

shinyUI(fluidPage(
  tags$head(
    tags$style(".resetButton {float:right; font-size:12px;}")
  ),

  titlePanel("Infection Pathways from Sequence Data"),

  sidebarLayout(

      sidebarPanel(
        uiOutput('patInput'),
        uiOutput('movInput'),
        uiOutput('aliInput'),
        uiOutput('wardInput'),
        HTML("<button id='reset' class='action-button resetButton'>Reset</button>")
      ),

      mainPanel(
          h4("Output here"),
          tableOutput(outputId = "reshead"),
          plotOutput(outputId = "rasterPlot")
      )
  )

))
