#****************************************************
#********      DASHBOARD: PGSMs for SBMs       ******
#****************************************************
rm(list = ls())
source("SourceScript.R"); set.seed(5)
#****************************************************

# generate SBM
num.nodes <- 20 # no. nodes
K <- 5  # no. clusters
pref.matrix <- diag(K) # Bernoulli rates (K x K matrix)
block.sizes <- rep(4, K) # no. nodes in each cluster (K length vector)
directed <- FALSE
sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
#start.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
#                       c(2, 6, 8, 10, 11, 17))
adj <- as_adj(sbm)

# PGSMs tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
N <- 20    # no. particles: (Bouchard uses 20)
resampling.threshold <- 0.5
n.iters <- 1000

# perform PGSMs
clusters <- start.clusters
num.clusters <- rep(0, n.iters)
start <- proc.time()
for(i in 1:n.iters)
{
  s <- SelectAnchors(clusters)
  clusters <- SplitMerge(s, clusters, adj, N, resampling.threshold, alpha, beta1, beta2, directed)
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
  num.clusters[i] <- length(clusters)
  print(num.clusters[i])
}
(run.time <- proc.time() - start)
