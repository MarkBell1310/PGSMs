#****************************************************
#********    Perform PGSMs for SBMs        **********
#****************************************************
rm(list = ls())
library(Matrix)
library(igraph)
library(matrixStats)
library(LaplacesDemon)
source("PGSMsFunctions.R")
#****************************************************

# generate SBM
n <- 20 # no. nodes
K <- 4  # no. clusters
directed <- FALSE
#pref.matrix <- diag(K) # Bernoulli rates (K x K matrix)
#block.sizes <- rep(4, K) # no. nodes in each cluster (K length vector)
#sbm <- sample_sbm(n, pref.matrix, block.sizes, directed = FALSE, loops = FALSE); plot(sbm)
set.seed(2)
sbm <- sample_sbm(n = 20,
                  pref.matrix = forceSymmetric(matrix(runif(K^2), (c(K, K)))),
                  block.sizes = c(6, 4, 7, 3), directed = FALSE, loops = FALSE); plot(sbm)
start.clusters <- list(1:n) # list of clusters
adj <- as_adj(sbm)

# PGSMs tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
N <- 5    # no. particles: (Bouchard uses 20)
resampling.threshold <- 0.5
n.iters <- 100000

# perform PGSMs
clusters <- start.clusters
num.clusters <- rep(0, n.iters)
start <- proc.time()
for(i in 1:n.iters)
{
  clusters <- SplitMerge(clusters, adj, N, resampling.threshold, alpha, beta1, beta2)
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
  num.clusters[i] <- length(clusters)
  print(num.clusters[i])
}
(run.time <- proc.time() - start)

