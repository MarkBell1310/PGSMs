#original.clusters <- all.clusters
#cluster.history <- lapply(1:n.iters, list)

all.clusters <- start.clusters
all.clusters

s <- SelectAnchors(all.clusters)  # select anchors uniform
s

# calculate c.bar and s.bar
closure <- CalculateClosureOfAnchors(s, all.clusters)
c.bar <- closure$c.bar
s.bar <- closure$s.bar

# calculate non.c.bar
non.c.bar.cluster.indicators <- lapply(lapply(all.clusters, 
                                              function(x){s %in% x}), 
                                       function(y){any(y == TRUE)})
non.c.bar <- all.clusters[which(non.c.bar.cluster.indicators == FALSE)]

# uniform permutation on elements of s.bar - with anchors 1st and 2nd
sigma <- SamplePermutation(s, s.bar)
n <- length(sigma)

# define particle matrix & fix 1st particle of each generation to conditional path
particles <- matrix(rep(0, N*n), c(N, n))
particles[1,] <- MapClustersToAllocations(sigma, c.bar) # row 1: particle with conditional path
particles[,1] <- rep(1, N) # column 1: 1st decision is (#1) initialise for all particles

# define weights at t=1
log.un.weights <- rep(log(1), N)
log.norm.weights <- rep(log(1/N), N) # 1st log norm weight is log(1/N) for all particles
log.weights <- matrix(rep(0, N*n), c(N, n)); log.weights[,1] <- log.un.weights


# Run SMC
for(t in 2:n) # time iterations
{
  # resample if ESS too small
  if(ESS(log.norm.weights, N) < resampling.threshold)
  {
    # resample and change history of resampled particles
    particles <- ResampleAndChangeParticleHistory(log.norm.weights, particles, N, t-1)
    log.un.weights <- rep(log(1), N) # reset the weights
  }
  
  for(p in 1:N) # particle iterations
  {
    # proposal allocation
    if(p >= 2)
    {
      # proposal only uses particle up to time t-1
      particles[p, t] <- Proposal(sigma, s, particles[p, 1:(t-1)], all.clusters, non.c.bar, 
                                  adj, tau1, tau2, t, n, alpha, beta1, beta2)
    }
    
    # unnormalised weights
    log.un.weights[p] <- LogUnnormalisedWeight(sigma, s, particles[p, 1:t], 
                                               log.previous.weight = log.un.weights[p], 
                                               all.clusters, non.c.bar, adj, tau1, tau2, 
                                               t, n, alpha, beta1, beta2)
  }
}

# normalised weights
log.norm.weights <- sapply(log.un.weights, function(x){x - logSumExp(log.un.weights)})

particles # 7 7
exp(log.norm.weights)

# choose a particle by multinomial sampling according to the weights
multinomial.sample <- as.double(rmultinom(n = 1, size = 1, prob = exp(log.norm.weights)))
updated.c.bar <- MapAllocationsToClusters(sigma, particles[which(multinomial.sample == 1), ], s)
updated.c.bar

unchanged.clusters <- all.clusters[all.clusters %!in% c.bar]

if(is.null(updated.c.bar[[2]]))
{
  # if merge we only have 1 updated cluster
  updated.clustering <- c(list(as.integer(updated.c.bar[[1]])), 
                          unchanged.clusters)
}
if(!is.null(updated.c.bar[[2]]))
{
  updated.clustering <- c(updated.c.bar, unchanged.clusters)
}
all.clusters <- updated.clustering
all.clusters

