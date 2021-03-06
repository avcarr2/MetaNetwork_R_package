% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WGCNAFunctions.R
\name{tidyModuleDataForOutput}
\alias{tidyModuleDataForOutput}
\title{Tidies the module eigenprotein data for easy access by ggplot2.}
\usage{
tidyModuleDataForOutput(dataInput, moduleColors, groupsFile)
}
\arguments{
\item{dataInput}{untidy data as output by the module eigenproteins workflow.}

\item{moduleColors}{merged module colors}

\item{groupsFile}{user-uploaded file with experimental groups and sampleIDs}
}
\value{
outputs "tidy" module data for easy plotting using ggplot2.
}
\description{
Tidies the module eigenprotein data for easy access by ggplot2.
}
