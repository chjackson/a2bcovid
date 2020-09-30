OurApp <- function(dat, rstudio=FALSE){
    launch.browser <- if (!rstudio) TRUE else rstudioapi::viewer
    shiny::runApp(system.file("OurApp", package = "a2bcovid"),
                  launch.browser = launch.browser)
}
