% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentWorkflowFunctions.R
\name{enrichAllez}
\alias{enrichAllez}
\title{Performs Allez enrichment}
\usage{
enrichAllez(
  GeneSymbols,
  GeneUniverse,
  SpeciesLibrary,
  idtype = "SYMBOL",
  alter = TRUE,
  Lowersetsize = 5,
  Uppsersetsize = 500,
  outprefix,
  ...
)
}
\arguments{
\item{GeneSymbols}{Gene symbols}

\item{GeneUniverse}{The total genes identified in this experiment}

\item{SpeciesLibrary}{The species corresponding to the samples.}

\item{idtype}{Defaults to "SYMBOL," meaning the standard gene symbol.}

\item{alter}{Defaults to true. Filters high and low set size values.}

\item{Lowersetsize}{Default 5. Terms below this count are filtered out .}

\item{Uppsersetsize}{Default 500. Terms above this count are filtered out.}

\item{outprefix}{The file written.}

\item{...}{Passes arguments}
}
\value{
Returns Allez gene ontology enrichment results.
}
\description{
Performs Allez enrichment
}
