#****************************************************
#************  DASHBOARD: Gibbs sampler  ************
#****************************************************
rm(list = ls())
source("SourceScript.R"); set.seed(5)
#****************************************************

# generate SBM
num.nodes <- 20 # no. nodes
#num.clust <- 4  # no. clusters
#pref.matrix = forceSymmetric(matrix(rbeta(num.clust^2, 2, 2),(c(num.clust, num.clust))))
#block.sizes <- rep(num.nodes/num.clust, num.clust) # no. nodes in each cluster (K length vector)
directed <- FALSE
sbm <- sample_sbm(n = 20, pref.matrix = diag(5),
                  block.sizes = rep(4, 5), directed = FALSE, loops = FALSE); plot(sbm)
#sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
adj <- as_adj(sbm)
all.clusters <- start.clusters

# Gibbs tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
n.iters <- 100000
prior <<- "mcdaid" # choose from "dirichlet.process" or "mcdaid"

# setup initialisation vectors/matrices
previous.matrices <- suppressMessages(InitialSetupList(all.clusters = start.clusters, 
                                                       num.nodes, adj, alpha, beta1, beta2, 
                                                       K = length(start.clusters), directed))

#****************************************************
# Run Gibbs sampler

clustering.list <- sapply(1:n.iters, function(x){list(x)})
total.cost <- 0
num.clusters <- rep(0, n.iters)
start <- proc.time()
for(i in 1:n.iters)
{
  run.Gibbs <- suppressMessages(GibbsSweepSolo(all.clusters, alpha, beta1, beta2, num.nodes, 
                                               previous.matrices, adj, directed = directed))
  all.clusters <- run.Gibbs$all.clusters
  previous.matrices <- run.Gibbs$previous.matrices
  
  # Record output  
  total.cost <- total.cost + (num.nodes * (length(all.clusters))^2)
  clustering.list[[i]] <- all.clusters
  num.clusters[i] <- length(all.clusters)
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
  print(num.clusters[i])
}
(run.time <- proc.time() - start) 
 
