
#****************************************************
#*************    Source script   *******************
#****************************************************

library(Matrix)
library(igraph)
library(matrixStats)
library(LaplacesDemon)
library(ergm)
library(Rcpp)
library(RcppArmadillo)
source("Functions - Initialisation.R")
sourceCpp("FunctionsCppGibbs.cpp")
source("Functions - PGSMs.R")
source("Functions - Gibbs.R")

#****************************************************