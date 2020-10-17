library(a2bcovid)
library(shiny)
library(shinyBS)
library(ggplot2)

defaultpars <- formals(a2bcovid) #list(seq_noise = 0.772469, evo_rate = 0.0008, min_qual = 0.8, max_n = 10)
