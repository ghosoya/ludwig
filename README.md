# ludwig
This repository contains the sourcecode of the R-package "ludwig".
The intended functionality of ludwig is to estimate the weights and
thresholds of an undirected network given a matrix of binary states.
Internally, ludwig uses the R-packages glmnet and glm.

In addition, it provides some functionality to visualize the network 
structure using qgraph as well as some means to simulate data from 
a network using conditional distributions.

New in Version 0.0.2

Bayesian estimation (MCMC) via JAGS as well as functionality for 
cross-validating estimated networks has been added.

Please note that this is a development version.
