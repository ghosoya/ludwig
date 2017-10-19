
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
#' @param X a binary data matrix containing observed network states. 
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
#'   \item {\code{id_index}} {a vectorindicating the states in \code{y}}.
#'   \item {\code{link_index}} {a  matrix that identifies the links between the nodes}.
#'   \item {\code{n_links}} {the number of links in the network}.
#'   \item {\code{n_param}} {the number of parameters to be estimated}.
#'   }
#'   
create_matrix<-function(X)
{
 size=dim(X) # dimensionality of matrix
 n=size[1] # number of state vectors
 k=size[2] # number of nodes
  
 n_links=0.5*(k^2-k) # number of links
 n_param=n_links+k # number of parameters
 res<-NULL
 
 # Initialization of dummy matrix for parameters
 interaction_matrix<-matrix(rep(0, n_links*n), nrow=n, ncol=n_links)
 # Initialization of link index
 link_index<-matrix(rep(NA, n_links*2), nrow=n_links, ncol=2)
 link_index
 dim(interaction_matrix)
 
 index<-1
 # Setting up the link index for each node
 for(i in 1:(k-1))
 {
   for(j in (i+1):k)
   {
     link_index[index,1]=i;
     link_index[index,2]=j;
 #### DEBUG cat(i, " ", j, "\n");
     index=index+1
   }
 }
 link_index
 
 # Constructing the predictor matrix
 predict_matrix<-c()
 counter<-0
## pb<-txtProgressBar(style=3, max=k*(k-1))
 for(i in 1:k) # k is number of nodes
 {
   
   node=i
   
   intercept_matrix<-matrix(rep(0, n*k), nrow=n, ncol=k)
   intercept_matrix[,node]=rep(1,n)
   interaction_matrix_copy<-interaction_matrix
  
   # cat("yay!")
   select1=which(link_index[,1]==node)
   parent1<-link_index[select1,2]
   select2=which(link_index[,2]==node)
   parent2<-link_index[select2,1]
   
   select<-c(select1, select2)
   parents<-c(parent1, parent2)
  
    for(j in 1:length(select))
   {
     interaction_matrix_copy[,select[j]]=X[,parents[j]]
     counter=counter+1
    # setTxtProgressBar(pb, counter,3)
   }
   
   predict_matrix_node<-cbind(intercept_matrix, interaction_matrix_copy)
   predict_matrix<-rbind(predict_matrix, predict_matrix_node)
 }
#### PB cat("\n")
 
 # Constructing the vector of observed values
 y<-c()
 
 for(i in 1:k)
 {
   y<-c(y,X[,i])
 }
 
 # Create ID-Index
 id_index<-c();
   for(i in 1:k)
   {
      id_index<-c(id_index, seq(1:n))
   }
 
 res$data<-X # Original data matrix
 res$predict_matrix<-predict_matrix # Design matrix
 res$link_index<-link_index # link index
 res$y<-y # observed values
 res$n<-n # number of observed state vectors
 res$k<-k # number of nodes
 res$id_index<-id_index# index of state vectors
 res$n_links<-n_links# links
 res$n_param<-n_param# number of parameters
 res
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
plot.estnet<-function(x, ...)
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
print.estnet<-function(x,...)
{
  cat("Thresholds:\n")
  cat(x$thresholds, "\n\n")
  cat("Weights:\n")
  print(x$weights)
  cat("\n Number of nodes: ", x$k)
  cat("\n Number of state vectors: ", x$n, "\n")
  cat("\n Number of data points:", length(x$y))
  cat("\n Number of parameters:", x$n_param)
  cat("\n log-Pseudolikelihood: ", x$logL)
  cat("\n Deviance: ", x$deviance)
  cat("\n")
  cat("\n Method:", x$method)
  if(x$method=="glmnet")
  {
  cat("\n Alpha: ", x$alpha)
  cat("\n s:", x$s)
  }
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
#'    (argument \code{s="lambda.min"} passed to function \code{\link[glmnet]{predict.cv.glmnet}}).
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
estnet<-function(data, method="glmnet",alpha=0.5, s="lambda.min",...)
{
  start.time <- Sys.time()
  res<-create_matrix(data)
  X<-res$predict_matrix
  y<-res$y
  
  if(method=="glmnet") # default: Use elastic-Net with cross validation
  {
    net = cv.glmnet(x=X, y=y, family="binomial", alpha=alpha,intercept=FALSE)
    coef<-coef(net, s=s)[1:res$n_param+1,1]
    res$method<-"glmnet"
    res$alpha<-alpha
    res$s<-s
  }
  else if(method=="glm") # use good old glm()
  {
    net<-glm(res$y~-1+res$predict_matrix, family="binomial")
    coef<-coef(net)
    res$method<-"glm"
    res$alpha<-NULL
    res$s<-NULL
  }
  else
  {
    cat("The argument you provided for method is not recognized.")
    return(-1)
  }
  # Link index
  link_index<-res$link_index
  # Create weight matrix
  weights<-matrix(NA, res$k*res$k, nrow=res$k, ncol=res$k)
  diag(weights)<-0
  for(i in 1:dim(res$link_index)[1])
  {
    weights[link_index[i,1], link_index[i,2]]=coef[i+res$k]
    weights[link_index[i,2], link_index[i,1]]=coef[i+res$k]
  }

  # Predicted probabilities 
  res$fitted<-exp(X%*%coef)/(1+exp((X%*%coef)))
  res$logL<-sum(log(exp(X%*%coef)^res$y/(1+exp((X%*%coef)))))
  #res$deviance<-deviance.glmnet(m1)[length(deviance.glmnet(m1))]
  res$deviance<--2*res$logL
  res$coefficients<-coef
  res$weights<-weights
  res$thresholds<-coef[1:res$k]
  class(res)<-"estnet"
  end.time <- Sys.time()
  res$time<-end.time-start.time
  res$net<-net
  res
}
#estnet<-function(data, method="glmnet", alpha=0.5, s="lambda.min", ...) UseMethod("estnet")
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
simnet<-function(k_sim, coef, n_burnin, n_rep, start_state, order="random")
{
  
  # initial start sate is random
  if(missing(start_state))
  {
    start_state=round(runif(k_sim))
  }
  
  # initial node configuration given by user
  data_sim<-matrix(start_state, nrow=1, ncol=k_sim) 
  create_matrix(data_sim)
  data_rep<-NULL

  # Predictor vector for initial node configuration
  res_sim<-create_matrix(data_sim)
  X<-res_sim$predict_matrix
  
  if(order=="sequential")
  {
    cat("Sequential updating\n")
    # Burn-in
    for(j in 1:n_burnin)  
    {
      for(i in 1:dim(X)[1])
      {
        lin_pred<-X[i,]%*%coef 
        p<-exp(lin_pred)/(1+exp(lin_pred))
        data_sim[i]=runif(1)<p
        res_sim<-create_matrix(data_sim)
        X<-res_sim$predict_matrix
      }
    }

    # Sample from the network
    for(j in 1:(n_rep))  
    {
      for(i in 1:dim(X)[1])
      {
        lin_pred<-X[i,]%*%coef 
        p<-exp(lin_pred)/(1+exp(lin_pred))
        data_sim[i]=runif(1)<p
        res_sim<-create_matrix(data_sim)
        X<-res_sim$predict_matrix
      }
      data_rep<-rbind(data_rep, data_sim)
    }
  data_rep
  }
  else if(order=="random")
  {
    cat("Random updating\n")
    for(j in 1:n_burnin)
    {
      #Choose one node at random for updating
      node<-sample(1:k_sim,1)
      lin_pred<-X[node,]%*%coef
      p<-exp(lin_pred)/(1+exp(lin_pred))
      data_sim[node]=runif(1)<p
      res_sim<-create_matrix(data_sim)
      X<-res_sim$predict_matrix
    }
    for(j in 1:n_rep)
    {
      #Choose one node at random for updating
      node<-sample(1:k_sim,1)
      lin_pred<-X[node,]%*%coef
      p<-exp(lin_pred)/(1+exp(lin_pred))
      data_sim[node]=runif(1)<p
      res_sim<-create_matrix(data_sim)
      X<-res_sim$predict_matrix
      data_rep<-rbind(data_rep, data_sim)
    }
    data_rep
  }
}
