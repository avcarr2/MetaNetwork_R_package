% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WGCNAFunctions.R
\name{mergeAndWriteWGCNAWorkbook}
\alias{mergeAndWriteWGCNAWorkbook}
\title{Merges modules broken up by list with the complete dataset sheet.}
\usage{
mergeAndWriteWGCNAWorkbook(selectedDatabase, allData, dataList, ...)
}
\arguments{
\item{selectedDatabase}{User uploaded database of UniProt Accessions and Genes}

\item{allData}{The result of the WGCNA}

\item{dataList}{List containing module and expression data}

\item{...}{passes arguments to workbook generation function}
}
\value{
This function combines allData, containing the complete set of module data,
with the dataList object, including a list with each index corresponding to a
different module.
}
\description{
Merges modules broken up by list with the complete dataset sheet.
}
