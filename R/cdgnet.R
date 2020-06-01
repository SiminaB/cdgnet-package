#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()

  cat("list paths within app function")
  print(names(list_paths_KEGG))
  shinyApp(ui = .cdgnetUI, server = .cdgnetServer)
}

#' run the CDGnet app
#'
#' @note Make sure to run function [`download_and_process_KEGG`] before running the app
#' @usage runApp()
#'
#' @export
runCDGnet <- function() {
  this_dir <- getwd()
  on.exit(setwd(this_dir))

  setwd(file.path(system.file("appdir", package="CDGnet")))
  source("app.R")
  shiny::runApp(cdgnetApp())
}
