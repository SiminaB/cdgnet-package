#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(list_paths_KEGG_file, database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()

  load(list_paths_KEGG_file)
  cat("list paths within app function")
  print(names(list_paths_KEGG))

#  e <- environment(.cdgnetServer)
#  assign("list_paths_KEGG", list_paths_KEGG, e)

#  stopifnot(exists("list_paths_KEGG", envir=e))

  print(ls(environment(.cdgnetServer)))
  shinyApp(ui = .cdgnetUI, server = eval(.cdgnetServer, envir=e))
}

#' run the CDGnet app
#'
#' @note Make sure to run function [`download_and_process_KEGG`] before running the app
#' @usage runApp()
#'
#' @export
runCDGnet <- function() {
  list_paths_KEGG_file <- file.path(system.file("appdir", package="CDGnet"), "list_paths_KEGG.RData")

  if (!file.exists(list_paths_KEGG_file)) {
    stop("KEGG data not found. Download with function CDGnet::download_and_process_KEGG()")
  }


  shiny::runApp(cdgnetApp(list_paths_KEGG_file))
}
