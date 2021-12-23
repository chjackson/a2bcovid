a2bcovid
===

The development repository for the `a2bcovid` R package.

## Installation (development version)

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/a2bcovid")
 ```

`a2bcovid` is a tool to estimate the likelihood that an infection was transmitted between particular individuals, and then to identify clusters of infections. 

It uses data on either genome sequences or the locations of infected individuals, and combines them in a statistical and evolutionary model.

In short our method compares data from a pair of individuals, A and B, to an underlying model of direct COVID transmission from A to B.  The philosophy of our method is based upon excluding cases: Rather than aiming to say whether A infected B, we evaluate whether, given the data, A *could have* infected B, excluding cases where the data are not consistent with direct transmission.

The key output from our method is a numerical value.  Given data from two individuals, A and B, we calculate the log probability of having observed these data given that A infected B.

The log probability is used to calculate a p-value for rejecting the null hypothesis of direct transmission.  If the p-value is greater than 0.05 we cannot reject this hypothesis, and say that the data is consistent with direct transmission from A to B (or in other words, A *could have* infected B).  If the p-value is less than 0.01, we reject the hypothesis, saying that direct transmission from A to B is unlikely.  If the p-value is intermediate between these numbers, we flag the case as being 'borderline'.

It is currently designed to apply to COVID-19 infection dynamics on hospital wards.
  
See the package vignette for an example:

[A2BCovid: Example analysis](https://chjackson.github.io/a2bcovid/docs/example.html)

Source code [GitHub repository](https://github.com/chjackson/a2bcovid)
