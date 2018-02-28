#****************************************************
#********    Perform PGSMs for SBMs        **********
#****************************************************
rm(list = ls())
library(Matrix)
library(igraph)
library(matrixStats)
source("PGSMsFunctions.R")
#****************************************************

# generate SBM
num.nodes <- 20
#num.clusters <- 4
# sbm <- sample_sbm(n = 20,
#                   pref.matrix = matrix(c(0.99, 0, 0, 0.1), c(2, 2)),
#                   block.sizes = c(10, 10), directed = FALSE, loops = FALSE); plot(sbm)
# sbm <- sample_sbm(n = 20,
#                   pref.matrix = forceSymmetric(matrix(runif(num.clusters^2),
#                                                       (c(num.clusters, num.clusters)))),
#                   block.sizes = c(6, 4, 7, 3), directed = FALSE, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20,
                  pref.matrix = forceSymmetric(
                    matrix(c(1, 0.01, 0.01, 1), c(2, 2))),
                  block.sizes = c(10, 10), directed = FALSE, loops = FALSE); plot(sbm)
# sbm <- sample_sbm(n = 20, pref.matrix = diag(5),
#                   block.sizes = rep(4, 5), directed = FALSE, loops = FALSE); plot(sbm)
# sbm <- sample_sbm(n = 20, pref.matrix = 0.99 * diag(1),
#                   block.sizes = 20, directed = FALSE, loops = FALSE); plot(sbm)

all.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
                     c(2, 6, 8, 10, 11, 17))
#all.clusters <- list(1:10, 11:num.nodes)
#all.clusters <- list(1:num.nodes)
adj <- as_adj(sbm)

# tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
N <- 20     # Bouchard also uses 20 particles
resampling.threshold <- 0.5
n.iters <- 100000 # (try 1000 and if run times are OK go higher)

# perform PGSMs
clusters <- all.clusters
num.clusters <- rep(0, n.iters)
start <- proc.time()
for(i in 1:n.iters)
{
  clusters <- SplitMerge(clusters, adj, N, resampling.threshold, alpha, beta1, beta2)
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
  num.clusters[i] <- length(clusters)
  print(num.clusters[i])
}
proc.time() - start
