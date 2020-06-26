#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(kegg_path, database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()

  shinyApp(ui = .cdgnetUI,
           server = .cdgnetServer,
           onStart = function() {
             assign("list_paths_KEGG",
                    readRDS(kegg_path),
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
  .check_KEGG(package=TRUE)
  kegg_path <- .keggFile(package=TRUE)
  shiny::runApp(cdgnetApp(kegg_path))
}
