% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{generateFixed}
\alias{generateFixed}
\title{Generate simulation predictors}
\usage{
generateFixed(P, start = 1)
}
\arguments{
\item{P}{number of predictors to simulate.}

\item{start}{numberic, where to begin the formula. By default, 1}
}
\value{
list of instructions to generate the predictors.
}
\description{
Indicates how to create each predictor to the simulation algorithm.
It will produce a list of predictors with the instructions to be generated
as a normal-like distribution with mean 0 and standard deviation 1.
}
