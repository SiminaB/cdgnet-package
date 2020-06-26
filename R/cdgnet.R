#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()
  .check_KEGG()

  shinyApp(ui = .cdgnetUI,
           server = .cdgnetServer,
           onStart = function() {
             assign("list_paths_KEGG",
                    readRDS(.keggFile()),
                    envir= globalenv() ) }
  )
}

#' run the CDGnet app
#'
#' @note Make sure to run function [`download_and_process_KEGG`] before running the app
#' @usage runApp()
#'
#' @export
runCDGnet <- function() {
  shiny::runApp(cdgnetApp())
}
