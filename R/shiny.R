a2bcovid_app <- function(dat, rstudio=FALSE){
    launch.browser <- if (!rstudio) TRUE else rstudioapi::viewer
    shiny::runApp(system.file("a2bcovid", package = "a2bcovid"),
                  launch.browser = launch.browser)
}
