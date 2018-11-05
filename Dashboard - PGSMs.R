#****************************************************
#**************** DASHBOARD: PGSMs ******************
#****************************************************
rm(list = ls())
source("SourceScript.R"); set.seed(5)
#****************************************************

## generate SBM
num.nodes <- 20 # no. nodes
num.clust <- 4  # no. clusters
pref.matrix = forceSymmetric(matrix(rbeta(num.clust^2, 2, 2),(c(num.clust, num.clust))))
block.sizes <- rep(num.nodes/num.clust, num.clust) # no. nodes in each cluster (K length vector)
directed <- FALSE
#sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20, pref.matrix = diag(5),
                  block.sizes = rep(4, 5), directed = FALSE, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
adj <- as_adj(sbm)
all.clusters <- start.clusters
global.num.clusters <<- length(all.clusters) # updated in SplitMerge()

## PGSMs tuning parameters
alpha <- 1 # Dirichlet process parameter
beta1 <- 1 # Flat uniform priors (McDaid: conjugate priors on the parameter for each cluster)
beta2 <- 1
N <- 20   # no. particles: (Bouchard uses 20)
resampling.threshold <- 0.5
n.iters <- 10000
as.probability <- 0.0 # probability of ancester sampling step 
prior <<- "dirichlet.process" # choose from "dirichlet.process" or "mcdaid"
network.model <<- "sbm"  # choose from "sbm" or "ldergm"
# model.formula <<- "edges + kstar(2)" # for LDERGM
# model.formula.terms <<- 2 # for LDERGM
# within.cluster.parameter.prior.density <<- dnorm # for LDERGM

## setup initialisation vectors/matrices
previous.matrices <- suppressMessages(InitialSetupList(all.clusters = start.clusters, 
                                                       num.nodes, adj, alpha, beta1, beta2, 
                                                       K = length(start.clusters), directed))
## perform PGSMs
num.clusters <- rep(0, n.iters)
clustering.list <- sapply(1:n.iters, function(x){list(x)})
total.cost <- 0
start <- proc.time()
for(i in 1:n.iters)
{
  PGSMs <- SplitMerge(all.clusters, adj, N, resampling.threshold, alpha, beta1, beta2, 
                      directed, as.probability, previous.matrices)
  all.clusters <- PGSMs$updated.clustering
  previous.matrices <- PGSMs$previous.matrices
  total.cost <- total.cost + PGSMs$num.nodes.in.c.bar * length(all.clusters) * N 
  
  # Record output  
  clustering.list[[i]] <- all.clusters
  global.num.clusters <<- num.clusters[i] <- length(all.clusters)
  print(num.clusters[i])
  if(i %% 10 == 0) {cat(paste0("iteration: ", i, "\n"))}
}
(run.time <- proc.time() - start)

