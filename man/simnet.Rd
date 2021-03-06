% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ludwig_functions.R
\name{simnet}
\alias{simnet}
\title{Simulate data from a network}
\usage{
simnet(k_sim, coef, n_burnin, n_rep, start_state, order = "random")
}
\arguments{
\item{k_sim}{the number of nodes.}

\item{coef}{a vector containing the coefficients of the network to simulate from. 
The first parameters in the vector are the thresholds. The remainig 
parameters are the upper diagonal of a weight matrix in sequential, row-wise order.}

\item{n_burnin}{the number of burn-in cycles.}

\item{n_rep}{the desired number of sampled states.}

\item{start_state}{provide a binary verctor as a starting state. If omitted, 
a random state will be generated as a starting state.}

\item{order}{either \code{random} (default) for randomly updating the network one node at a time or \code{sequential} 
for sequential updating. With \code{random} updating, each new state will be returned. 
If \code{sequential} is used, a state is returned after all \code{k_sim} nodes have been
updated.}
}
\description{
This function is intended to simulate data from a binary network given
    a user defined or randomly generated starting state and a vector of coefficients.
    Two methods of simulating states from the conditional distributions are provided: 
    random updating or sequential updating.
    A burn-in parameter is provided to allow for 
    running the update process a desired number of iterations before the 
    actual sampling process starts.
}
\examples{
# Simulates two negatively coupled nodes
# using random updating
simnet(2, c(1,1,-2), 1000, 1000)
}
