#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()
  shinyApp(ui = .cdgnetUI, server = .cdgnetServer)
}

#' run the CDGnet app
#'
#' @usage runApp()
#'
#' @export
runApp <- function() {
  source(file.path(system.file("appdir", package="CDGnet"), "app.R"))
}
