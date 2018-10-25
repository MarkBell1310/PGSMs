#****************************************************
#*****  DASHBOARD: PGSMs with Gibbs iterations ******
#****************************************************
rm(list = ls())
source("SourceScript.R"); set.seed(5)
#****************************************************

## Generate SBM
num.nodes <- 20 # no. nodes
K <- 4  # no. clusters
pref.matrix = forceSymmetric(matrix(rbeta(K^2, 2, 2),(c(K, K))))
block.sizes <- rep(num.nodes/K, K) # no. nodes in each cluster (K length vector)
directed <- FALSE
sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
adj <- as_adj(sbm)
all.clusters <- start.clusters

## Tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
N <- 20   # no. particles: (Bouchard uses 20)
resampling.threshold <- 0.5
n.iters <- 100000
as.probability <- 0.0 # probability of ancester sampling step 
gibbs.frequency <- 1 # frequency of performing Gibbs sweep
prior <<- "mcdaid" # choose from "dirichlet.process" or "mcdaid"
network.model <<- "sbm"  # choose from "sbm" or "ldergm"
# model.formula <<- "edges + kstar(2)" # for LDERGM
# model.formula.terms <<- 2 # for LDERGM
# within.cluster.parameter.prior.density <<- dnorm # for LDERGM

## Perform PSGMs/Gibbs split merge algorithm
num.clusters <- rep(0, n.iters)
clustering.list <- sapply(1:n.iters, function(x){list(x)})
total.cost <- 0
start <- proc.time()
for(i in 1:n.iters)
{
  # PGSMs
  if(i %% 2 == 1)
  {
    PGSMs <- SplitMerge(all.clusters, adj, N, resampling.threshold, alpha, beta1, beta2, 
                        directed, as.probability)
    all.clusters <- PGSMs$updated.clustering
    total.cost <- total.cost + PGSMs$num.nodes.in.c.bar * length(all.clusters) * N 
    
    # Record output  
    clustering.list[[i]] <- all.clusters
    num.clusters[i] <- length(all.clusters)
    if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
    print(num.clusters[i])
  }

  # Gibbs sampler
  #if(i %% gibbs.frequency == 0)
  if(i %% 2 == 0)
  {
    all.clusters <- suppressMessages(GibbsSweepIntegrated(all.clusters, alpha, beta1, beta2, num.nodes, 
                                                          adj, directed = directed))
    total.cost <- total.cost + num.nodes * (length(all.clusters))^2 
    
    # Record output  
    clustering.list[[i]] <- all.clusters
    num.clusters[i] <- length(all.clusters)
    if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
    print(num.clusters[i])
  }
}
(run.time <- proc.time() - start)


