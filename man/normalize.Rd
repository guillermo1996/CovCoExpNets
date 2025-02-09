% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{normalize}
\alias{normalize}
\title{Normalize a vector}
\usage{
normalize(x)
}
\arguments{
\item{x}{numeric vector or list of numeric vectors}
}
\value{
numeric vector or list of numeric vectors
}
\description{
Normalizes a vector by subtracting its mean and dividing by its standard
deviation
}
\examples{
normalize(1:4)
}
