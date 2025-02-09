% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{removeSexualGenes}
\alias{removeSexualGenes}
\title{Remove sexual genes}
\usage{
removeSexualGenes(data)
}
\arguments{
\item{data}{numerical matrix with predictors as rows and samples as columns.}
}
\value{
reduced numerical matrix
}
\description{
Identifies and eliminates genes located in the sexual and mithocondrial chromosomes.
}
