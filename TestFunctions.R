
#****************************************************
#
#**** Particle Gibbs for SBMs - Test functions ******
#
#****************************************************
rm(list = ls())
library(Matrix)
library(igraph)
library(matrixStats)
#library(blockmodels)
source("PGSMsFunctions.R")
#****************************************************

# generate SBM
num.nodes <- 20
num.clusters <- 4
sbm <- sample_sbm(n = num.nodes, 
                  pref.matrix = forceSymmetric(matrix(runif(num.clusters^2), 
                                                      (c(num.clusters, num.clusters)))), 
                  #pref.matrix = matrix(c(0.99, 0.01, 0.01, 0.99), c(num.clusters, num.clusters)),
                  block.sizes = c(6, 4, 7, 3),
                  directed = FALSE,
                  loops = FALSE)
plot(sbm)
adj <- as_adj(sbm)

# define clusters: 20 data points in 4 clusters
all.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15), 
                     c(2, 6, 8, 10, 11, 17))

# select anchors at random
s <- SelectAnchors(all.clusters)   

# calculate c.bar and s.bar
closure <- CalculateClosureOfAnchors(s, all.clusters)
c.bar <- closure$c.bar
s.bar <- closure$s.bar

# uniform permutation on closure of anchors
sigma <- SamplePermutation(s, s.bar)         

# get particle from c.bar
particle <- MapClustersToAllocations(sigma, c.bar)
MapAllocationsToClusters(sigma, particle, s)        # should have same elements as c.bar
c.bar

# intermediate targets
alpha <- 1
beta1 <- 0.1
beta2 <- 0.2
t <- 3
n <- length(particle)
log.previous.weight <- log(0.1)

LogIntermediateTarget(sigma, s, particle, all.clusters, c.bar, adj, tau1, tau2, 
                      t, alpha, beta1, beta2)

LogImprovedIntermediateTarget(sigma, s, particle, all.clusters, c.bar, adj, tau1, tau2,
                              t, n, alpha, beta1, beta2)

# proposal and weights
PossibleAllocations(sigma, s, particle, all.clusters, c.bar, adj, 
                    tau1, tau2, t, n, alpha, beta1, beta2)

Proposal(sigma, s, particle, all.clusters, c.bar, adj, tau1, tau2,
         t, n, alpha, beta1, beta2)

LogUnnormalisedWeight(sigma, s, particle, log.previous.weight, all.clusters, 
                      c.bar, adj, tau1, tau2, t, n, alpha, beta1, beta2)

# PGSM
N <- 10
resampling.threshold <- 0.5

ParticleGibbsSplitMerge(all.clusters, adj, s, s.bar, c.bar, N, resampling.threshold,
                        alpha, beta1, beta2)

SplitMerge(all.clusters, adj, N, resampling.threshold, alpha, beta1, beta2)















