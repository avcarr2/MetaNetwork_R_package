## Call Shiny App ##
##' Starts MetaNetwork
#' @name RunMetaNetwork
#' @description Function to run MetaNetwork shiny application
#' @return Runs the MetaNetwork shiny application
#' @param ... Arguments to pass to shiny::runApp
#' @export

RunMetaNetwork <- function(results_directory){
  withr::with_dir(results_directory,
                    shiny::runApp(system.file("App",
                                              package = "MetaNetwork"))
                  )
}

