###########################################################
#' Create the predictor matrix and the vector of 
#' observed node states for pseudolikelihood estimation
#' 
#' @description The main purpose of this function is to create 
#'     a predictor matrix (\code{predict_matrix}) and a vector of observed states 
#'     (\code{y}) for each node
#'     so that the network parameters are estimable using simple
#'     multivariate logistic regression based on the conditional distributions
#'     of each node given the rest of the network. This method is also 
#'     known as the pseudolikelihood method, supposedly first introduced by Besag (1975).
#'     Note that the speed of this function is improvable.
#' @param data a binary data matrix containing observed network states. 
#'     The rows represent the states, the columns represent the
#'     nodes.
#' 
#' @return A list containing... 
#'   \itemize{
#'   \item {\code{data}} {the original data matrix}.
#'   \item {\code{predict_matrix}} {the predictor matrix}.
#'   \item {\code{y}} {a vector of the observed node states}.
#'   \item {\code{n}} {the number of observed state vectors}.
#'   \item {\code{k}} {the number of nodes}.
#'   \item {\code{n_links}} {the number of links in the network}.
#'   \item {\code{n_param}} {the number of parameters to be estimated}.
#'   }
#'   
create_matrix <- function(data)
{
  data_size <- length(data)
  
  n <- dim(data)[1L]
  k <- dim(data)[2L]
  
  n_links <- as.integer(k * (k - 1L) / 2L)
  
  predict_matrix <- matrix(0L, nrow = data_size, ncol = n_links + k)

  for (j in seq_len(k))
  {
    predict_matrix[(n * (j - 1L) + 1L) : (n * j), j] <- 1L
  }
  
  col_idx = k + 1L
  
  for (j_1 in seq_len(k - 1L))
  {
    for (j_2 in (j_1 + 1L) : k)
    {
      predict_matrix[(n * (j_1 - 1L) + 1L) : (n * j_1), col_idx] <- data[, j_2]
      predict_matrix[(n * (j_2 - 1L) + 1L) : (n * j_2), col_idx] <- data[, j_1]
      
      col_idx = col_idx + 1L
    }
  }
  
  res_list <- list(data = data,
                   predict_matrix = predict_matrix,
                   y = c(data),
                   n = n,
                   k = k,
                   n_links = n_links,
                   n_param = n_links + k)
  
  return(res_list)
}
###########################################################
#' Plot the estimated network
#' 
#' @description This function plots the estimated network.
#'     Basically, it is just a simple wrapper around  \code{\link[qgraph]{qgraph}},
#'     so all of qgraph's beautiful functionality should be available.
#' 
#' @param x an object of class \code{estnet}
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{\link[qgraph]{qgraph}}
#' @export
plot.estnet <- function(x, ...)
{
  qgraph(x$weights, ...)
}
############################################################
#' Print the results 
#' 
#' @description This function prints the weights and
#'     thresholds of the estimated network as well as 
#'     other useful information.
#'     
#' @param x An object of class \code{estnet}
#' @param ... other arguments (not used)
#' @export
print.estnet <- function(x, ...)
{
  cat("Thresholds:\n")
  cat(x$thresholds, "\n\n")
  cat("Weights:\n")
  print(x$weights)
  cat("\n\n")
  cat("\n Number of nodes: ", x$k)
  cat("\n Number of state vectors: ", x$n, "\n")
  cat("\n Number of data points:", length(x$y))
  cat("\n Number of parameters:", x$n_param)
  cat("\n log-Pseudolikelihood: ", x$logL)
  cat("\n Deviance: ", x$deviance)
  cat("\n")
  cat("\n Method:", x$method)
  if (x$method == "glmnet")
  {
    cat("\n Alpha: ", x$alpha)
    cat("\n s:", x$s)
  }
  cat("\n\n Entropy rate (bits):", x$entropy_rate)
  cat("\n\n Time:", x$time)
}
##########################################################
#' Estimate the network parameters given an observed state matrix 
#'
#'@description This function accepts a binary matrix of network states and estimates
#'    the thresholds and weights of the network model.
#'    The internally used functions are either \code{\link[glmnet]{cv.glmnet}} or \code{\link[stats]{glm}}.
#'    By default, \code{\link[glmnet]{cv.glmnet}} is used. 
#'    In this case, the mixing parameter \code{alpha} is set to \code{0.5} 
#'    (elastic net regularization) and \code{lambda} 
#'    is chosen based on the solution with mimimum mean cross-validation error 
#'    (argument \code{s = "lambda.min"} passed to function \code{\link[glmnet]{predict.cv.glmnet}}).
#'
#'@param data a binary matrix of observed network states. The rows are the states, the columns are the nodes.
#'@param method the method to be used. Either \code{"glmnet"} (default) or \code{"glm"}. 
#'@param alpha the elastic mixing parameter passed to \code{\link[glmnet]{cv.glmnet}}. 
#'@param s value of the penalty parameter passed to \code{\link[glmnet]{predict.cv.glmnet}}.
#'@param ... other arguments passed over to e.g. \code{\link[glmnet]{cv.glmnet}}.
#'@return The function returns a list of class \code{\link{estnet}}, containing...
#'   \itemize{
#'   \item {the objects returned by the function \code{\link{create_matrix}}}.
#'   \item {\code{fitted}} { fitted values based on the estimated parameters.}
#'   \item {\code{y}} { a vector of the observed node states.}
#'   \item {\code{logL}} { the log pseudo likelihood.}
#'   \item {\code{deviance}} { the deviance.}
#'   \item {\code{coefficients}} { the estimated coefficients (thresholds and weights).}
#'   \item {\code{weights}} { the estimated weights.}
#'   \item {\code{thresholds}} { the estimated thresholds.}
#'   \item {\code{time}} { the time used for the estimation process.}
#'   \item {\code{net}} { the list returned by either \code{\link[stats]{glm}} or \code{\link[glmnet]{cv.glmnet}}.}
#'   }
#'   
#'@seealso  \code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}, \code{\link[glmnet]{predict.cv.glmnet}},
#'     \code{\link[stats]{glm}}
#'@export
estnet <- function(data, method = "glmnet", alpha = 0.5, s = "lambda.min", ...)
{
  if (class(data) != "matrix")
  {
    stop("Input must be a matrix.")
  }
  
  if (!all(unique(c(data)) %in% c(0,1)))
  {
    stop("If the input matrix is not of class logical, it may only contain ones and zeros.")
  }
  
  storage.mode(data) <- "integer"
  
  start.time <- Sys.time()
  res <- create_matrix(data)

  if (method == "glmnet") # default: Use elastic-Net with cross validation
  {
    net <- cv.glmnet(x = res$predict_matrix, y = res$y, family = "binomial", alpha = alpha, intercept = FALSE)
    coef <- coef(net, s = s)[seq_len(res$n_param) + 1L, 1L]
    res$method <- "glmnet"
    res$alpha <- alpha
    res$s <- s
  }
  else if (method == "glm") # use good old glm()
  {
    net <- glm(res$y ~ -1 + res$predict_matrix, family = "binomial")
    coef <- coef(net)
    res$method <- "glm"
    res$alpha <- NULL
    res$s <- NULL
  }
  else
  {
    stop("The argument you provided for method is not recognized.")
  }
  # Create weight matrix
  weights <- matrix(as.numeric(0), nrow = res$k, ncol = res$k)
  coef_idx <- res$k + 1L
  
  for (j_1 in seq_len(res$k - 1L))
  {
    for (j_2 in (j_1 + 1L) : res$k)
    {
      weights[j_1, j_2] <- coef[coef_idx]
      weights[j_2, j_1] <- coef[coef_idx]
      
      coef_idx <- coef_idx + 1L
    }
  }
  # Predicted probabilities 
  res$fitted <- exp(res$predict_matrix %*% coef) / (1 + exp(res$predict_matrix %*% coef))
  res$entropy_rate<-mean(-((res$fitted)*log(res$fitted,2)+(1-res$fitted)*log(1-res$fitted,2)))
  res$logL <- sum(log(exp(res$predict_matrix %*% coef) ^ res$y / (1 + exp(res$predict_matrix %*% coef))))
  #res$deviance <- deviance.glmnet(m1)[length(deviance.glmnet(m1))]
  res$deviance <- -2 * res$logL
  res$coefficients <- coef
  res$weights <- weights
  res$thresholds <- coef[seq_len(res$k)]
  res$time <- Sys.time() - start.time
  res$net <- net
  
  class(res) <- "estnet"
  
  return(res)
}
##########################################################
#' Simulate data from a network 
#' 
#' @description This function is intended to simulate data from a binary network given
#'     a user defined or randomly generated starting state and a vector of coefficients.
#'     Two methods of simulating states from the conditional distributions are provided: 
#'     random updating or sequential updating.
#'     A burn-in parameter is provided to allow for 
#'     running the update process a desired number of iterations before the 
#'     actual sampling process starts.
#'     
#' @param k_sim the number of nodes.
#' @param coef a vector containing the coefficients of the network to simulate from. 
#'     The first parameters in the vector are the thresholds. The remainig 
#'     parameters are the upper diagonal of a weight matrix in sequential, row-wise order.  
#' @param n_burnin the number of burn-in cycles.
#' @param n_rep the desired number of sampled states.
#' @param start_state provide a binary verctor as a starting state. If omitted, 
#'    a random state will be generated as a starting state.
#' @param order either \code{random} (default) for randomly updating the network one node at a time or \code{sequential} 
#'     for sequential updating. With \code{random} updating, each new state will be returned. 
#'     If \code{sequential} is used, a state is returned after all \code{k_sim} nodes have been
#'     updated.
#' 
#' @examples
#' # Simulates two negatively coupled nodes
#' # using random updating
#' simnet(2, c(1,1,-2), 1000, 1000)
#' @export
simnet <- function(k_sim, coef, n_burnin, n_rep, start_state, order = "random")
{
  
  # initial start sate is random
  if (missing(start_state))
  {
    start_state <- round(runif(k_sim))
  }
  
  # initial node configuration given by user
  data_sim <- matrix(start_state, nrow = 1L)
  data_rep <- NULL
  
  # Predictor vector for initial node configuration
  res_sim <- create_matrix(data_sim)

  if (order == "sequential")
  {
    cat("Sequential updating\n")
    # Burn-in
    for (j in seq_len(n_burnin))
    {
      for (i in seq_len(res_sim$k))
      {
        lin_pred <- res_sim$predict_matrix[i, ] %*% coef
        p <- exp(lin_pred) / (1 + exp(lin_pred))
        data_sim[i] <- runif(1L) < p
        res_sim <- create_matrix(data_sim)
      }
    }
    # Sample from the network
    for (j in seq_len(n_rep))  
    {
      for (i in seq_len(res_sim$k))
      {
        lin_pred <- res_sim$predict_matrix[i, ] %*% coef 
        p <- exp(lin_pred) / (1 + exp(lin_pred))
        data_sim[i] <- runif(1L) < p
        res_sim <- create_matrix(data_sim)
      }
      data_rep <- rbind(data_rep, data_sim)
    }
    return(data_rep)
  }
  else if (order == "random")
  {
    cat("Random updating\n")
    # Burn-in
    for (j in seq_len(n_burnin))
    {
      #Choose one node at random for updating
      node <- sample(seq_len(k_sim), 1)
      lin_pred <- res_sim$predict_matrix[node, ] %*% coef
      p <- exp(lin_pred) / (1 + exp(lin_pred))
      data_sim[node] <- runif(1L) < p
      res_sim <- create_matrix(data_sim)
    }
    for (j in seq_len(n_rep))
    {
      #Choose one node at random for updating
      node <- sample(seq_len(k_sim), 1)
      lin_pred <- res_sim$predict_matrix[node, ] %*% coef
      p <- exp(lin_pred) / (1 + exp(lin_pred))
      data_sim[node] <- runif(1L) < p
      res_sim <- create_matrix(data_sim)
      data_rep <- rbind(data_rep, data_sim)
    }
    return(data_rep)
  }
}

##########################################################
#' Cross-validate network on new data
#' 
#' @description This intended use of this function is helping cross-validate
#'     an estimated network on new data. To this aim, the values of individual
#'     nodes are predicted based on the conditional predictive distributions
#'     and the relative frequency of correct predictions is returned.
#'
#' @param net An object of class estnet
#' @param dat_val A validation dataset
#' 
#' @return The relative frequency of correct predictions.

crossval <- function(net, dat_val)
{
  res <- create_matrix(dat_val)
  coef <- coef(net)
  pred <- round(exp(res$predict_matrix %*% coef) / (1 + exp(res$predict_matrix %*% coef)))
  return(mean(pred == res$y))
}

##########################################################
#' Determine the posterior distributions of the network parameters given an observed state matrix 
#' 
#' @description This function allows for estimating an undirected graphical model
#' using MCMC via JAGS (Just another Gibb's sampler). 
#'     
#' @param data a binary input matrix
#' @param n_chains the number of chains used
#' @param n_iter the number of MCMC iterations
#' @param n_burnin the number of burn-in iterations
#' @param n_adapt the number of adaption iterations
#' @param n_thin the thinning interval
#'   
#' @return The function returns a list of class \code{\link{estnet_bayes}}, containing...
#'   \itemize{
#'   \item {the objects returned by the function \code{\link{create_matrix}}}
#'   \item {\code{s_coda}} {the raw coda samples}
#'   \item {\code{fitted}} { expected values of the nodes given the posterior means of the network parameters}
#'   \item {\code{entropy_rate}} {the entropy rate based on the network parameters' posterior means}
#'   \item {\code{logL}} { the log pseudo likelihood}
#'   \item {\code{deviance}} { the deviance}
#'   \item {\code{merged chains}} { the merged MCMC chains}
#'   \item {\code{used_samples}} { the actual number of samples to assess the posterior distributions}
#'   \item {\code{posterior_means}} { the parameter's posterior means}
#'   \item {\code{posterior_sd}} { the parameter's posterior standard deviations}
#'   \item {\code{weights}} { a matrix of the network parameter's posterior means}
#'   \item {\code{weights_sd}} { a matrix of the network parameter's posterior standard deviations}
#'   \item {\code{weights}} { the estimated weights}
#'   \item {\code{thresholds}} { posterior means of the threshold parameters}
#'   \item {\code{thresholds_sd}} { posterior standard deviations of the threshold parameters}
#'   \item {\code{n_chains}} {number of MCMC chains}
#'   \item {\code{thinning}} {thinning interval}
#'   \item {\code{burnin}} {burnin iterations}
#'   \item {\code{n_iter}} {number of iterations}
#'   \item {\code{actual}} {actual iterations used to compute the posterior statistics}
#'   \item {\code{time}} { the time used for the estimation process}
#'   }
#'   
#' @examples
#' # Simulates two negatively coupled nodes
#' # using random updating
#' simnet(2, c(1,1,-2), 1000, 1000)
#' @export
estnet_bayes <- function(data, n_chains = 2, n_iter = 1000, n_burnin= 1000, n_adapt = 1000, n_thin = 4)
{
  if (class(data) != "matrix")
  {
    stop("Input must be a matrix.")
  }

  if (!all(unique(c(data)) %in% c(0,1)))
  {
    stop("If the input matrix is not of class logical, it may only contain ones and zeros.")
  }
  
  storage.mode(data) <- "integer"
  load.module("glm")
  
  start.time <- Sys.time()
  # Create design matrix
  res <- create_matrix(data)
  
  # Prepare data for jags
  y <- res$y
  X <- res$predict_matrix
  n <- length(y)
  p <- ncol(X) 
  
  data <- list("y" = y, "X" = X, "n" = n, "p" = p)
  
  # Parameters of interest
  parameters <- c("beta")
  
  # Inits drawn from a standard normal distribution
  inits <- function() { list(beta = rnorm(p)) }
  
  jags <- jags.model(file = textConnection(log_reg),
                     data = data,
                     inits = inits,
                     n.chains = n_chains,
                     n.adapt = n_adapt)
  
  update(jags, n_adapt)
  
  # Coda Samples
  s_coda <- coda.samples(jags,
                         parameters,
                         n.iter = n_iter,
                         thin = n_thin)
  
  # Merge chains
  merged_chains <- mcmc(do.call(rbind, s_coda))
  # Compute posterior M
  posterior_means <- apply(merged_chains, 2, mean)
  
  # Create weight matrix
  weights <- matrix(as.numeric(0), nrow = res$k, ncol = res$k)
  coef_idx <- res$k + 1L
  
  for (j_1 in seq_len(res$k - 1L))
  {
    for (j_2 in (j_1 + 1L) : res$k)
    {
      weights[j_1, j_2] <- posterior_means[coef_idx]
      weights[j_2, j_1] <- posterior_means[coef_idx]
      
      coef_idx <- coef_idx + 1L
    }
  }
  
  # Compute posterior SD
  posterior_sd <- apply(merged_chains, 2, sd)
  
  # Posterior posterior SD matrix
  weights_sd <- matrix(as.numeric(0), nrow = res$k, ncol = res$k)
  coef_idx <- res$k + 1L
  
  for (j_1 in seq_len(res$k - 1L))
  {
    for (j_2 in (j_1 + 1L) : res$k)
    {
      weights_sd[j_1, j_2] <- posterior_sd[coef_idx]
      weights_sd[j_2, j_1] <- posterior_sd[coef_idx]
      
      coef_idx <- coef_idx + 1L
    }
  }
  
  res$s_coda <- s_coda # Raw coda output
  res$fitted <- exp(res$predict_matrix %*% posterior_means) / (1 + exp(res$predict_matrix %*% posterior_means))
  res$entropy_rate<-mean(-((res$fitted)*log(res$fitted,2)+(1-res$fitted)*log(1-res$fitted,2)))
  res$logL <- sum(log(exp(res$predict_matrix %*% posterior_means) ^ res$y / (1 + exp(res$predict_matrix %*% posterior_means))))
  res$deviance <- -2 * res$logL
  res$merged_chains <- mcmc(do.call(rbind, s_coda)) # Combined chains
  res$used_samples <- dim(merged_chains) # Actual Samples
  res$posterior_means <- apply(merged_chains, 2, mean)
  res$posterior_sd <- apply(merged_chains, 2, sd)
  res$weights <- weights # weight posterior means matrix
  res$weights_sd <- weights_sd # weights posterior sd matrix
  res$thresholds <- posterior_means[seq_len(res$k)] # Thresholds posterior means
  res$thresholds_sd <- posterior_sd[seq_len(res$k)] # Thresholds posterior SD
  res$n_chains<-n_chains
  res$thinning<-n_thin
  res$burnin<-n_burnin
  res$iter<-n_iter
  res$actual<-dim(res$merged_chains)[1]
  res$time <- Sys.time() - start.time
  
  class(res) <- "estnet_bayes"
  return(res)
}

############################################################
#' Print the results 
#' 
#' @description This function prints the weights and
#'     thresholds of the estimated network as well as 
#'     other useful information.
#'     
#' @param x An object of class \code{estnet_bayes}
#' @param ... other arguments (not used)
#' @export
print.estnet_bayes <- function(x, ...)
{
  cat("Thresholds (Posterior Means):\n")
  cat(x$thresholds, "\n\n")
  cat("Thresholds (Posterior SD):\n")
  cat(x$thresholds_sd, "\n\n")
  cat("Weights (Posterior Means):\n")
  print(x$weights)
  cat("\n")
  cat("Weights (Posterior SD):\n")
  print(x$weights_sd)
  cat("\n\n")
  cat("\n Number of nodes: ", x$k)
  cat("\n Number of state vectors: ", x$n, "\n")
  cat("\n Number of data points:", length(x$y))
  cat("\n Number of parameters:", x$n_param)
  cat("\n log-Pseudolikelihood: ", x$logL)
  cat("\n Deviance: ", x$deviance)
  cat("\n")
  cat("\n Method: Bayes" )
  cat("\n Iterations:", x$iter )
  cat("\n Chains:", x$n_chains )
  cat("\n Burn in:", x$burnin )
  cat("\n Thinning:", x$thinning  )
  cat("\n Number of samples used:", x$actual  )
  cat("\n\n Entropy rate (bits):", x$entropy_rate)
  cat("\n\n Time:", x$time)
}
###########################################################
#' Plot the estimated network
#' 
#' @description This function plots the estimated network.
#'     Basically, it is just a simple wrapper around  \code{\link[qgraph]{qgraph}},
#'     so all of qgraph's functionality should be available.
#' 
#' @param x an object of class \code{estnet_bayes}
#' @param thresh a threshold in terms of posterior mean divided by posterior sd. 
#'     For instance, setting the threshold to 2 removes all edges for which the 
#'     absolute qoutient of posterior mean and posterior sd is smaller than 2.
#'     This is useful for removing edges with small parameter values.
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{\link[qgraph]{qgraph}}
#' @export
plot.estnet_bayes <- function(x, thresh = 0, ...)
{
  qgraph(x$weights, ...)
  cat("Bayes")
  
  if(thresh!=0)
  {
      w <- x$weights*(x$weights/x$weights_sd>threshold|x$weights/x$weights_sd< -threshold)
      diag(w)<-0
      qgraph(w, ...)
  }
}
##########################################################
#' Determine the convergence of the Markov chains
#' 
#' @description This function is a wrapper around \code{gelman.diag} to be found 
#'     in the \code{coda} package. It computes the potential
#'     scale reduction for the Markov chains. This statistic is needed for the 
#'     assessment of model convergence.
#' 
#' @param x an object of class \code{estnet_bayes}
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{gelman.diag}, \code{coda}
#' @export
R_hat<-function(x, ...)
{
  gelman.diag(x$s_coda, ...)
}

##########################################################
#' Check autocorrelations
#'
#' @description This function is a wrapper around \code{autocorr.diag} to be found 
#'     in the \code{coda} package. It computes the autocorrelations
#'     for the Markov chains
#' 
#' @param x an object of class \code{estnet_bayes}
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{autocorr.diag}, \code{coda}
#' @export
auto<-function(x, ...)
{
 autocorr.diag(x$s_coda, ...) 
}
#########################################################
#' Diagnostic plots
#'
#' @description This function is a wrapper around \code{mcmcplots} to be found 
#'     in the \code{coda} package. It allows for a convenient visual 
#'     assessment of the Markov chains and the autocorrelation in a Web browser.
#' 
#' @param x an object of class \code{estnet_bayes}
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{autocorr.diag}, \code{coda}
#' @export
diag_plot<-function(x, ...)
{
  mcmcplot(x$s_coda, ...)
}
#########################################################
#' Deviance information criterion
#'
#' @description This function is a wrapper around \code{DIC} to be found 
#'     in the \code{coda} package. It allows for assessing the 
#'     deviance information criterion.
#' 
#' @param x an object of class \code{estnet_bayes}
#' @param ... other arguments passed over to \code{qgraph}
#' 
#' @seealso \code{autocorr.diag}, \code{coda}
#' @export
DIC<-function(x, ...)
{
  dic.samples(x$s_coda, ...)
}

########################################################
#' Basic logistic regression model
#' 

log_reg<-"model{
  
  # Likelihood  
  for(i in 1:n)
  {
    logit(prob[i])<-inprod(X[i,], beta)
    y[i]~dbern(prob[i])
  }
  
  # Priors
  for(j in 1:p)
  {
    beta[j]~dnorm(0, 0.5)
  }
}"