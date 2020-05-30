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
