% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ludwig_functions.R
\name{auto}
\alias{auto}
\title{Check autocorrelations}
\usage{
auto(x, ...)
}
\arguments{
\item{x}{an object of class \code{estnet_bayes}}

\item{...}{other arguments passed over to \code{autocorr.diag()}}
}
\description{
This function is a wrapper around \code{autocorr.diag} to be found 
    in the \code{coda} package. It computes the autocorrelations
    for the Markov chains
}
\seealso{
\code{autocorr.diag}, \code{coda}
}
