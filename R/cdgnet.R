## Class constructor
cdgnetApp <- function(database=NULL, names=NULL) {
  # Run the application
  shinyApp(ui = .cdgnetUI, server = .cdgnetServer)
}
