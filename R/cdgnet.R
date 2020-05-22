#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(database=NULL, names=NULL) {
  .check_and_load_KEGG()

  # Run the application
  shinyjs::useShinyjs()
  shinyApp(ui = .cdgnetUI, server = .cdgnetServer)
}
