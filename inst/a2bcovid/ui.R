txtfiles_accepted <-   c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                         'text/plain', '.csv', '.tsv')

shinyUI(fluidPage(
  titlePanel("Infection Pathways from Sequence Data"),
  
  sidebarLayout(
    
      sidebarPanel(
      fileInput('pat',
                'Patient file in CSV format. [ Info about format here ].',
                accept = txtfiles_accepted),
      fileInput('mov',
                'Movement file in CSV format. [ Info about format here ].',
                accept = txtfiles_accepted),
      fileInput('ward',
                'Ward file in CSV format. [ Info about format here ].',
                accept = txtfiles_accepted),
      fileInput('ali',
                'Sequence file in FASTA format. [ Info about format here ].',
                accept = ".fa")
      ),
      
      mainPanel(
          h4("Output here"),
          tableOutput(outputId = "reshead"),
          plotOutput(outputId = "rasterPlot")
      )
  )

))
