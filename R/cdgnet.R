#' Class constructor
#' @import shiny
#'
#'@export
#'
cdgnetApp <- function(list_paths_KEGG_file, database=NULL, names=NULL) {
  # Run the application
  shinyjs::useShinyjs()

  data(DrugBank_targets)
  data(drugs_PO_FDA_biomarkers)
  data(drugs_PO_FDA_targets)
  data(KEGG_cancer_paths_onc_long)
  data(FDA_approved_drugs)
  data(non_prop2prop)
  data(prop2non_prop)
  data(prop_non_prop)
  data(Onc_df)

  load(list_paths_KEGG_file)
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
  list_paths_KEGG_file <- file.path(system.file("appdir", package="CDGnet"), "list_paths_KEGG.RData")

  if (!file.exists(list_paths_KEGG_file)) {
    stop("KEGG data not found. Download with function CDGnet::download_and_process_KEGG()")
  }


  shiny::runApp(cdgnetApp(list_paths_KEGG_file))
}
