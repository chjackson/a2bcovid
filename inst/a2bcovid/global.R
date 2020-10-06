library(a2bcovid)
library(shiny)
library(ggplot2)

data_choices <- c("Just times of symptom onset",
  "Times of symptom onset and genome sequence data",
  "Times of symptom onset, sequence data, and patient locations",
  "Times of symptom onset, sequence data, patient locations and staff locations")
