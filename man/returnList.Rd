% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{returnList}
\alias{returnList}
\title{Return function for all list functions}
\usage{
returnList(return.list, data)
}
\arguments{
\item{return.list}{boolean. Whether to return a list or not}

\item{data}{list}
}
\value{
the list or the first element of it
}
\description{
Interal function to return either the whole list or only the first element.
It is employed when a function can return both a single item or a list of items
(usually, when several conditions are processed at the same time)
}
