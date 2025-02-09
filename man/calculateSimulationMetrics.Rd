% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{calculateSimulationMetrics}
\alias{calculateSimulationMetrics}
\title{Calculates different simulation model's metrics}
\usage{
calculateSimulationMetrics(set.calculated, sim.obj, covariate = "Q1")
}
\arguments{
\item{set.calculated}{vector of predictors selected by the model.}

\item{sim.obj}{simulation object obtained from CovCoExpNets.}

\item{covariate}{which covariate to study. Either "Q1" or "Q2". Defaults
to Q2.}
}
\value{
simulation model's metrics.
}
\description{
Given a set of predictors selected by the model (set.calculated)
and the simulation object to extract the total number of predictors,
this function calculates the precision, sensitivity, specificity, jaccard
index and the returned predictors.
}
