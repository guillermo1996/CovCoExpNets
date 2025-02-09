% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{plotRMSE}
\alias{plotRMSE}
\title{Plots the RMSE}
\usage{
plotRMSE(
  df.group,
  output.file = NA,
  title = "",
  ylab = "",
  legend.title = "Algorithm",
  first.arg = "GLMNET",
  second.arg = "CovCoExpNets"
)
}
\arguments{
\item{df.group}{data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}}

\item{output.file}{path the output file to save the plot. If not provided, the plot will
not be stored}

\item{title}{string, title of the plot}

\item{ylab}{string, y-axis title}

\item{legend.title}{string, title of the legend}
}
\value{
ggplot object containing the RMSE comparison of the two groups
}
\description{
For the given grouped results (see \link[CovCoExpNets]{groupMetrics}),
it plots the RMSE for the conditions and its 95\% confidence
interval.
}
