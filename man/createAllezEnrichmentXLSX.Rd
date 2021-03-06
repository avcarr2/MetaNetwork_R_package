% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentWorkflowFunctions.R
\name{createAllezEnrichmentXLSX}
\alias{createAllezEnrichmentXLSX}
\title{Writes the Allez Enrichment excel worksheet to the results folder}
\usage{
createAllezEnrichmentXLSX(geneUniverse, AllezPvalues)
}
\arguments{
\item{geneUniverse}{A vector containing all gene symbols used as the background in enrichment}

\item{AllezPvalues}{A list containing the calculated p-values from the addPvaluesToAllezOutput function}
}
\value{
Writes the results of the Allez enrichment workbook in .xlsx format
}
\description{
Writes the Allez Enrichment excel worksheet to the results folder
}
