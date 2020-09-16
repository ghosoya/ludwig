#' : Simulating and Estimating Undirected Binary Networks
#' 
#'     
#'     The purpose of this package is to provide some tools for methodological
#'     researchers interested in exploring undirected binary networks.
#'     The parameter estimation method is inspired by an article by Strauss (1992).
#'     In this article, Strauss describes with reference to Besag (1975), how lattice models could be fit using
#'     the pseudolikelihood method. As a consequence, binary network models could be fit
#'     with standard statistical functions capable of multivariate logistic regression, 
#'     such as \code{glmnet} or \code{glm}. 
#'     This package provides an experimental infastructure for exploring this idea
#'     with possible extensions in mind. However, it has to 
#'     be noted that the sampling properties of the maximum pseudolikelhood 
#'     estimator (MPE) seem to be unexplored (see Strauss, 1992). 
#'     Similar to the package \code{IsingFit} by van Borkulo, Epskamp and Robitzsch (2014), this package uses 
#'     regularized logistic regression (package \code{\link[glmnet]{glmnet}}) by default. In addition, the use of 
#'     the standard package \code{\link[stats]{glm}} is available. 
#'     
#'     Note that Bayesian estimation is only available when 
#'     JAGS \href{http://mcmc-jags.sourceforge.net.}{http://mcmc-jags.sourceforge.net.} is installed.
#' 
#' 
#' @section Functions:
#' 
#'   \itemize{
#'    \item {\code{\link{estnet}}} {Estimate the network parameters using regularized logistic regression using \code{\link[glmnet]{glmnet}}}.
#'    \item {\code{\link{simnet}}} {Simulate data from a network using probabilistic sequential node updating}.
#'    \item {\code{\link{create_matrix}}} {An internal function to create the predictor matrix and data vector for maximum pseudolikelihood estimation}.
#'    \item {\code{\link{print.estnet}}} {Print the results of the analysis}.
#'    \item {\code{\link{plot.estnet}}} {Plot the results of the analysis using \code{\link[qgraph]{qgraph}}}.
#'    \item {\code{\link{crossval}}} {A crossvalidation function to assess the network's predictive capabilities on the global and node levels}.
#'    \item {\code{\link{estnet_bayes}}} {Bayesian estimation of the network parameters using \code{rjags}}.
#'    \item {\code{\link{print.estnet_bayes}}} {Prints the results of the Bayesian estimation}.
#'    \item {\code{\link{plot.estnet_bayes}}} {Plot the results of the analysis using \code{\link[qgraph]{qgraph}}}.
#'    \item {\code{\link{R_hat}}} {Prints the Gelman-Rubin-Statistic for the Markov chains via \code{gelman.diag}}.
#'    \item {\code{\link{auto}}} {Prints autocorrelation of the chains via \code{autocorr.diag}}.
#'    \item {\code{\link{diag_plot}}} {Visualizes the chains and convergence diagnostics via \code{mcmcplot}}.
#'    \item {\code{\link{DIC}}} {Assesses the deviance information criterion (DIC) via \code{dic.samples}}.
#'   }
#' 
#' @section References:
#'
#' Besag, J. (1975). Statistical analysis of non-lattice data. \emph{The Statistician}, 24(\emph{3}), 179--195.
#'
#' Strauss, D. (1992). The many faces of logistic regression. \emph{American Statistician}, 46(\emph{4}), 321--327. 
#' 
#' van Borkulo, C., Epskamp, S., & Robitzsch (2014). IsingFit: Fitting Ising models using the eLasso method. R package version 0.3.0. 
#' 
#' @examples 
#' # Create the coefficients that define a
#' # network with four nodes
#' coef1 <- c(-0.5,-0.5,-0.5,-0.5,1,-1,1,1,-1,1)
#' 
#' # Simulate 1000 observations from the network.
#' # Use a burnin period of 5000 iterations and
#' # random updating.
#' dat1 <- simnet(4,coef1,5000,1000)
#'
#' # Try to recover the parameters
#' # using the default settings of estnet()
#' net1 <- estnet(dat1)
#'   
#' # Print the results
#' print(net1)
#'
#' # Try the same with glm()
#' net2 <- estnet(dat1, method="glm")
#'
#' # Print the results
#' print(net2)
#'
#' # Plot the networks
#' plot(net1, labels = c("A", "B", "C", "D"), 
#'      maximum = 1.5)
#' plot(net2, labels = c("A", "B", "C", "D"),
#'      maximum = 1.5) 
#'      
#' # Print estimated vs true coefficients
#' plot(coef(net1), coef1, xlab = "Estimated", ylab="True")
#' plot(coef(net2), coef1, xlab = "Estimated", ylab="True")
#'
#'
#' @docType package
#' @name ludwig
#' @import glmnet
#' @import qgraph
#' @import stats
#' @import mcmcplots
NULL