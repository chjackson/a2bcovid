##' Web app interface to a2bcovid
##'
##' @param rstudio Set to \code{TRUE} to open the app in the RStudio Viewer.
##' If \code{FALSE} (the default), an external web browser is launched.
##'
##' @export
a2bcovid_app <- function(rstudio=FALSE){
    launch.browser <- if (!rstudio) TRUE else rstudioapi::viewer
    shiny::runApp(system.file("a2bcovid", package = "a2bcovid"),
                  launch.browser = launch.browser)
}
