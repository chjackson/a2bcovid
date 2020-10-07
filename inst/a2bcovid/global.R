library(a2bcovid)
library(shiny)
library(ggplot2)

data_choices <- c("Just times of symptom onset" = "times",
  "Times of symptom onset and genome sequence data" = "times_gen",
  "Times of symptom onset, sequence data, and patient locations" = "times_gen_pat",
  "Times of symptom onset, sequence data, patient locations and staff locations" = "times_gen_pat_staff")

include_seq <- 2:4
include_patloc <- 3:4
include_staffloc <- 4

## could do
## box for include staff locations - requires patient locations and genome data
## box for include patient locations - requires genome data

# what to do with file dialogs then
# only include file dialog if
