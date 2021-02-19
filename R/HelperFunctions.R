## Helper Functions ##

#' Reads an excel sheet and imports each sheet as a component of a list.
#' @param filename path to excel sheet to be loaded
#' @param tibble make the data imported a tibble. Default = FALSE. Data imported as a data.frame.
#' @return List where each list elements corresponds to an excel sheet.
#' @name read_excel_allsheets
#' @importFrom utils read.csv
#' @export
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#' Override for read.csv that sets stringsAsFactors to TRUE.
#' @param fileName default file name to pass to results folder.
#' @param ... Params passed to read.csv. See documentation for read.csv for full list of parameters
#' @return override for read.csv that sets stringsAsFactors = TRUE by default. Will prevent future errors
#' @name read.csv_override
#' @export
read.csv_override <- function(fileName, ...){
  withr::with_options(
      list(stringsAsFactors = TRUE),
      read.csv(file = fileName, ...)
    )
}


# May need this function in the future
#  Uses withr to temporarily set the working directory:
#
# write_csv_toUserSelectedFile <- function(objectToWrite,
#                                          pathToWrite,
#                                          ...){
#   withr::with_dir(
#     pathToWrite,
#     write.csv(x = objectToWrite, file = pathToWrite)
#   )
# }
