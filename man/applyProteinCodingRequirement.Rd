% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{applyProteinCodingRequirement}
\alias{applyProteinCodingRequirement}
\title{Remove non protein-coding genes}
\usage{
applyProteinCodingRequirement(data, includeSexChromosomes = T)
}
\arguments{
\item{data}{numerical matrix with predictors as rows and samples as columns}

\item{includeSexChromosomes}{whether to include genes from sexual chromosomes
or not. By default, true.}
}
\value{
reduced numerical matrix
}
\description{
Identifies and eliminates genes which are not protein coding. It can also
remove genes from sexual chromosomes.
}
