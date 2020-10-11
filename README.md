a2bcovid
===

The development repository for the [a2bcovid](http://cran.r-project.org/package=fic) R package for 

## Installation (development version)

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/a2bcovid")
 ```

`a2bcovid` is a tool to estimate the likelihood that an infection was transmitted between particular individuals, and then to identify clusters of infections. 

It uses data on either genome sequences or the locations of infected individuals, and combines them in a statistical and evolutionary model.

It is currently designed to apply to COVID-19 infection dynamics on hospital wards.
  
See the package vignette for an example:

[A2BCovid: Example analysis](https://chjackson.github.io/a2bcovid/inst/doc/example.html)

Source code [GitHub repository](https://github.com/chjackson/a2bcovid)
