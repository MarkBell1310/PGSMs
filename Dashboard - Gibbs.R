 
#****************************************************
#************  DASHBOARD: Gibbs sampler *************
#****************************************************
rm(list = ls())
source("SourceScript.R"); set.seed(5)

#****************************************************
# generate SBM
num.nodes <- 20 # no. nodes
#K <- 4  # no. clusters
#pref.matrix <- diag(K) # Bernoulli rates (K x K matrix)
#pref.matrix = forceSymmetric(matrix(runif(K^2),(c(K, K))))
#block.sizes <- rep(num.nodes/K, K) # no. nodes in each cluster (K length vector)
directed <- FALSE
#sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
#start.clusters <- list(1:25, 26:50)
#start.clusters <- list(1:10, 11:20)
#start.clusters <- list(1:4, 5:8, 9:12, 13:16, 17:20)
#start.clusters <- list(1:10, 11:20, 21:30, 31:40, 41:50)
#start.clusters <- start.clusters <- list(1:5, 6:10, 11:15, 16:20, 21:25, 26:30, 31:35, 36:40, 41:45, 46:50)
#start.clusters <- sapply(1:num.nodes, function(x){list(x)})
# start.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
#                        c(2, 6, 8, 10, 11, 17))
sbm <- sample_sbm(n = num.nodes, pref.matrix = diag(5),
                 block.sizes = rep(4, 5), directed = directed, loops = FALSE); plot(sbm)
#sbm <- sample_sbm(n = num.nodes, pref.matrix = diag(2),
#                  block.sizes = rep(25, 2), directed = FALSE, loops = FALSE); plot(sbm)
adj <- as_adj(sbm)

# Gibbs tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
n.iters <- 100000

# setup initialisation vectors/matrices
previous.matrices <- suppressMessages(InitialSetupList(all.clusters = start.clusters, 
                                                       num.nodes, adj, alpha, beta1, beta2, 
                                                       K = length(start.clusters), directed))
all.clusters <- start.clusters

#****************************************************
# Run Gibbs sampler

num.clusters <- rep(0, n.iters)
start <- proc.time()
for(i in 1:n.iters)
{
  run.Gibbs <- suppressMessages(GibbsSweep(all.clusters, alpha, beta1, beta2, num.nodes, 
                                previous.matrices, directed = directed))
  all.clusters <- run.Gibbs$all.clusters
  previous.matrices <- run.Gibbs$previous.matrices
  
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
  num.clusters[i] <- length(all.clusters)
  print(num.clusters[i])
}
(run.time <- proc.time() - start) 
 
