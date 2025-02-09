% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{measureStabilityJaccard}
\alias{measureStabilityJaccard}
\title{Calculate similarity between glmnet repetitions}
\usage{
measureStabilityJaccard(genes.freq, diag = "default")
}
\arguments{
\item{genes.freq}{data.frame containing the GLMNET repetitions results.
Should be the output of \link[CovCoExpNets]{geneFrequency}.}

\item{diag}{(optinal) which parameter to register in the diagional. Only
other possible value is "genes". By default, the jaccard index.}
}
\value{
matrix with the similarity for between all glmnet repetitions
}
\description{
For a given \link[CovCoExpNets]{geneFrequency} object, it calculates
the jaccard index between all glmnet repetitions as a measure of
their similarity. Higher values represent more similar models and vice
versa.
}
