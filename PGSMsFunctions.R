
#****************************************************
#
#******* Particle Gibbs for SBMs - functions ********
#
#****************************************************

#****************************************************
#'  "Not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

#****************************************************
#'  Annealing schedule ("zeta") for improved target distributions at time t
#'  
#' @param t current time [scalar]
#' @param n maximum time [scalar]
#' @return annealing schedule [vector]
zeta_t <- function(t, n)
{
  (t - 2) / (n - 2)
}


#****************************************************
#'  Dirichlet process prior functions
tau1 <- function(alpha, j)
{
  alpha^j
}

tau2 <- function(j)
{
  factorial(j-1)
}


#****************************************************
#'  Effective sample size (ESS)
#'  @param log.norm.weights log of the normalised weights [vector]
#'  @param N number of particles [scalar]
#'  @return ESS
ESS <- function(log.norm.weights, N)
{
  ess <- exp(-log(N) - logSumExp(2*log.norm.weights))
  ifelse(is.nan(ess), 0, ess)
}

#****************************************************
#'  Select anchors
#'  
#' @param all.clusters The set of all clusters ("c") [list]
#' @return Two anchor points
SelectAnchors <- function(all.clusters)
{
  sample(as.double(unlist(all.clusters)), 2)
}


#****************************************************
#' Find Closure of the anchors
#'  
#' @param all.clusters The set of all clusters ("c") [list]
#' @param s Anchor points [vector]
#' @return c.bar [list] and s.bar [vector]  
CalculateClosureOfAnchors <- function(s, all.clusters)
{
  # find only the clusters that contain the anchors
  c.bar.index <- sapply(1:2, 
                        function(z)
                        {
                          as.double(which(sapply(all.clusters, function(x){s[z] %in% x}) == TRUE))
                        })
  
  c.bar <- list(all.clusters[[c.bar.index[1]]],
                all.clusters[[c.bar.index[2]]])
  
  # if both anchors in same cluster
  if(c.bar.index[1] == c.bar.index[2])
  {
    c.bar[[2]] <- NULL
  }
  
  # find all individual elements of c.bar
  s.bar <- as.double(unlist(c.bar))
  
  return(list("c.bar" = c.bar,
              "s.bar" = s.bar))
}
# CalculateClosureOfAnchors(s = c(2,1), all.clusters)


#****************************************************
#'  Perform uniform permutation on closure of anchors
#'  
#' @param s Anchors  [vector]
#' @param s.bar Closure of anchors [vector]
#' @return Uniform permutation [vector]  
SamplePermutation <- function(s, s.bar)
{
  # if 3 points in s.bar, then only 1 "non-anchor" point in s.bar
  if(length(s.bar == 3))
  {
    return(c(sample(s), s.bar[! s.bar %in% s]))
  }
  else
  {
    return(c(sample(s), sample(s.bar[! s.bar %in% s])))
  }
}
# # Test
# s <- c(6, 9)
# s.bar <- c(4, 7, 8, 1, 2, 6, 3, 9)
# SamplePermutation(s, s.bar)


#****************************************************
#'  Mapping from allocation decisions to clustering
#'  
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param s anchor points  [vector]
#' @return "c.bar": corresponding clustering [list]
MapAllocationsToClusters <- function(sigma, particle, s)
{
  # At t=2: if step 2 is merge then we merge throughout and return output
  if(particle[2] == 2)   
  {
    return(list(sigma, NULL))
  }
  
  # At t=2: otherwise split - assign an anchor to each cluster
  n <- length(sigma)
  cluster1 <- cluster2 <- rep(0, n) 
  cluster1[1] <- s[1]  # put each anchor back in the correct cluster
  cluster2[1] <- s[2]
  
  #cluster1[1] <- sigma[1] 
  #cluster2[1] <- sigma[2]
  
  # Then for t=3,...,n: calculate clustering from decisions #3 and #4
  if(n >= 3)
  {
    counter.c1 <- counter.c2 <- 2
    for(t in 3:n)
    {
      if(particle[t] == 3) # decision #3 adds to cluster 1
      {
        cluster1[counter.c1] <- sigma[t]
        counter.c1 <- counter.c1 + 1
      }
      if(particle[t] == 4) # decision #4 adds to cluster 2
      {
        cluster2[counter.c2] <- sigma[t]
        counter.c2 <- counter.c2 + 1
      }
    }
  }
  return(list(cluster1[cluster1!=0], cluster2[cluster2!=0]))
}
# ## Test
# p1 <- c(1, 4, 4, 3, 4, 3, 4, 4)
# p2 <- c(1, 4, 3, 3, 4, 3, 4, 3)
# p3 <- c(1, 2, 2, 2, 2, 2, 2, 2)
# MapAllocationsToClusters(sigma = s.bar, particle = p1)


#****************************************************
#'  Mapping from clustering to allocation decisions
#'  
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param c.bar The 1 or 2 clusters that contain the anchors [list]
#' @return "particle": corresponding allocation decison sequence [vector]
MapClustersToAllocations <- function(sigma, c.bar)
{
  n <- length(sigma)
  
  # if c.bar only has 1 cluster, anchors are in same cluster and hence
  # the allocation decision was merge
  if(length(c.bar)==1)
  {
    return(c(1, rep(2, n-1)))
  }
  
  # otherwise work backwards from clustering
  c1 <- c.bar[[1]]; len1 <- length(c1)
  c2 <- c.bar[[2]]; len2 <- length(c2)
  particle <- rep(0, n); particle[1:2] <- c(1, 4)
  
  for(t in 3:n)
  {
    # if data point is in cluster 1, assign allocation decision #3
    if(sigma[t] %in% c1) 
    {
      particle[t] <- 3
    }
    
    # if data point is in cluster 2, assign allocation decision #4
    if(sigma[t] %in% c2) 
    {
      particle[t] <- 4
    }
  }
  return(particle)
}
# ## Test
# c.bar <- MapAllocationsToClusters(sigma = s.bar, particle = p1)
# MapClustersToAllocations(sigma = s.bar, c.bar)
# p1


#****************************************************
#'  Count edges WITHIN the c.bar clusters
#'  (Assumes graph is undirected - and hence adjacency matrix is symmetric)
#'  
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @return list of "counts" and "max.counts" of edges  [list of vectors]
CountEdgesWithinCbarClusters <- function(c.bar.current, adj)
{
  # (allows for cases with 1 or 2 clusters in c.bar)
  counts.within.c.bar.clusters <- sapply(1:length(c.bar.current), function(x)
  {
    # multiply by 0.5 since we don't want to double count
    0.5 * sum(adj[c.bar.current[[x]], c.bar.current[[x]]]) 
  })
  
  max.counts.within.c.bar.clusters <- sapply(1:length(c.bar.current), function(x)
  {
    0.5 * length(c.bar.current[[x]]) * (length(c.bar.current[[x]]) - 1)
  })
  
  return(list("counts" = counts.within.c.bar.clusters,
              "max.counts" = max.counts.within.c.bar.clusters))
}

#****************************************************
#'  Count edges BETWEEN the c.bar clusters
#'  (Assumes graph is undirected - and hence adjacency matrix is symmetric)
#'  
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @return list of "counts" and "max.counts" of edges  [list of vectors]
CountEdgesBetweenCbarClusters <- function(c.bar.current, adj)
{
  # c.bar.current must have 2 clusters
  counts.between.c.bar.clusters <- sum(adj[c.bar.current[[1]], c.bar.current[[2]]])
  max.counts.between.c.bar.clusters <- length(c.bar.current[[1]]) * length(c.bar.current[[2]])
  
  return(list("counts" = counts.between.c.bar.clusters,
              "max.counts" = max.counts.between.c.bar.clusters))
}
  
#****************************************************
#'  Count edges BETWEEN the c.bar and non.c.bar clusters
#'  (Assumes graph is undirected - and hence adjacency matrix is symmetric)
#'  
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @return list of "counts" and "max.counts" of edges  [list of vectors]
CountEdgesBetweenCbarAndNonCbarClusters <- function(c.bar.current, adj, non.c.bar)
{
  counts.between.c.bar.non.c.bar.clusters <- as.double(sapply(non.c.bar, function(x)
  {
    if(length(c.bar.current) == 1) # if c.bar.current only contains 1 cluster
    {
      return(sum(adj[c.bar.current[[1]], x]))
    }
    if(length(c.bar.current) == 2)
    {
      return(c(sum(adj[c.bar.current[[1]], x]), sum(adj[c.bar.current[[2]], x])))
    }
  }))
  
  # max edge counts BETWEEN all c.bar.current and non.c.bar clusters
  max.counts.between.c.bar.non.c.bar.clusters <- as.double(sapply(non.c.bar, function(x)
  {
    if(length(c.bar.current) == 1) # if c.bar.current only contains 1 cluster
    {
      return(length(c.bar.current[[1]]) * length(x))
    }
    if(length(c.bar.current) == 2)
    {
      return(c(length(c.bar.current[[1]]) * length(x), length(c.bar.current[[2]]) * length(x)))
    }
  }))
  
  return(list("counts" = counts.between.c.bar.non.c.bar.clusters,
              "max.counts" = max.counts.between.c.bar.non.c.bar.clusters))
}
  
  
#****************************************************
#'  Count ALL relevant edges - used in the likelihood calculation
#'  (Assumes graph is undirected - and hence adjacency matrix is symmetric)
#'  
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @return list of "counts" and "max.counts" of edges  [list of vectors]
CountEdges <- function(c.bar.current, adj, non.c.bar)  
{
  # This function needs to consider the following scenarios:
  # (1) nodes WITHIN clusters in c.bar.current
  # (2) nodes BETWEEN clusters in c.bar.current (both clusters belong to c.bar) 
  # (3) nodes BETWEEN c.bar clusters and non.c.bar clusters

  # SCEN (1): nodes within clusters in c.bar
  edges.within.c.bar.current <- CountEdgesWithinCbarClusters(c.bar.current, adj)
  counts.within.c.bar.current <- edges.within.c.bar.current$counts
  max.counts.within.c.bar.current <- edges.within.c.bar.current$max.counts
  
  # SCEN (2): nodes between clusters in c.bar - c.bar must have 2 clusters
  if(length(c.bar.current) == 2)
  {
    edges.between.c.bar.current <- CountEdgesBetweenCbarClusters(c.bar.current, adj)
    counts.between.c.bar.current <- edges.between.c.bar.current$counts
    max.counts.between.c.bar.current <- edges.between.c.bar.current$max.counts
  }
  else
  {
    counts.between.c.bar.current <- NULL
    max.counts.between.c.bar.current <- NULL
  }
  
  # SCEN (3): nodes between c.bar and non c.bar clusters
  if(length(non.c.bar) != 0)
  {
    edges.between.c.bar.non.c.bar <- CountEdgesBetweenCbarAndNonCbarClusters(c.bar.current, 
                                                                             adj, non.c.bar)
    counts.between.c.bar.non.c.bar <- edges.between.c.bar.non.c.bar$counts
    max.counts.between.c.bar.non.c.bar <- edges.between.c.bar.non.c.bar$max.counts
  }
  else
  {
    counts.between.c.bar.non.c.bar <- NULL
    max.counts.between.c.bar.non.c.bar <- NULL
  }

  return(list("counts" = c(counts.within.c.bar.current,
                           counts.between.c.bar.current,
                           counts.between.c.bar.non.c.bar),
              "max.counts" = c(max.counts.within.c.bar.current,
                               max.counts.between.c.bar.current,
                               max.counts.between.c.bar.non.c.bar)))
}
#CountEdges(c.bar.current, adj, non.c.bar)  


#****************************************************
#'  Log of intermediate target distribution at time t (Log "Gamma")
#'
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param all.clusters set of all clusters ("c") [list]
#' @param c.bar.current clusters that contain the anchors, filled up to time t [list] 
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time
#' @param alpha tau1 parameter
#' @param beta1 beta function parameter
#' @param beta2 beta function parameter
#' @return log of intermeduate target distribution
LogIntermediateTarget <- function(sigma, s, particle, all.clusters, c.bar.current, non.c.bar, 
                                  adj, tau1, tau2, t, alpha, beta1, beta2)
{
  # count relevant edges within and between clusters
  edge.counts <- CountEdges(c.bar.current, adj, non.c.bar)  
  counts <- edge.counts$counts
  max.counts <- edge.counts$max.counts
  n <- length(counts)
  
  # product of all between-cluster log likelihoods 
  log.likelihoods <- rep(0, n)
  for (i in 1:n)
  {
    log.likelihoods[i] <-  log(beta(beta1 + counts[i], max.counts[i] - counts[i] + beta2)) - 
                            log(beta(beta1, beta2))
  }
  
  sum.log.likelihoods <- sum(log.likelihoods)
  
  # prior (3rd line is 0 if we split)
  log.prior <- log(tau1(alpha, length(c.bar.current) + length(all.clusters) - length(c.bar.current))) +
               log(tau2(length(c.bar.current[[1]]))) + 
               ifelse(particle[t] != 2, log(tau2(length(c.bar.current[[2]]))), 0) 

  # intermediate target
  log.int.target <- log.prior + sum.log.likelihoods
  
  return(log.int.target)
}   
#LogIntermediateTarget(sigma, s, particle, all.clusters, c.bar.current, non.c.bar, 
#                      adj, tau1, tau2, t, alpha, beta1, beta2)


#****************************************************
#' Log of Improved intermediate target distribution (Log "Gamma hat")
#'
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param all.clusters set of all clusters ("c") [list]
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time [scalar]
#' @param n maximum time [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter [scalar]
#' @param beta2 beta function parameter [scalar]
#' @return log of improved intermeduate target distribution
LogImprovedIntermediateTarget <- function(sigma, s, particle, all.clusters, non.c.bar, 
                                          adj, tau1, tau2, t, n, alpha, beta1, beta2)
{
  # if t = 1,2 then "gamma hat" = 1 and log "gamma hat" = 0
  if(t <= 2)
  {
    return(0)
  }
  
  # calculate "c.bar.current": c.bar at time t
  c.bar.current <- MapAllocationsToClusters(sigma[1:t], particle, s)
  if(is.null(c.bar.current[[2]]))
  {
    c.bar.current <- list(c.bar.current[[1]])
  }

  # log intermediate targets at current time t and time = 2
  log.gamma_t <- LogIntermediateTarget(sigma, s, particle, all.clusters, c.bar.current, non.c.bar, 
                                       adj, tau1, tau2, t, alpha, beta1, beta2)
  
  log.gamma_2 <- LogIntermediateTarget(sigma, s, particle, all.clusters, c.bar.current, non.c.bar, 
                                       adj, tau1, tau2, t = 2, alpha, beta1, beta2)
  
  # improved intermediate target
  log.improved.int.target <- (zeta_t(t, n) - 1) * log.gamma_2 + log.gamma_t
  
  return(log.improved.int.target)
}

#****************************************************
#'  Possible allocations - output used to calculate proposal and weights 
#'  -gives possible allocations - and the log gamma.hats of staying in current state or moving
#'  -function only called for t > 2
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t-1 [vector]
#' @param all.clusters set of all clusters ("c") [list]
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time [scalar]
#' @param n maximum time [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter [scalar]
#' @param beta2 beta function parameter [scalar]
#' @return the possible allocations and resulting log "gamma hats" at time t
PossibleAllocations <- function(sigma, s, particle, all.clusters, non.c.bar, 
                                adj, tau1, tau2, t, n, alpha, beta1, beta2)
{ 
  previous.allocation <- particle[t-1]
  
  # particle and log gamma.hat resulting from staying in current state 
  # (i.e. repeating previous allocation decision)
  stay.particle <- c(particle[1:(t-1)], previous.allocation)
  log.stay.gamma.hat <- LogImprovedIntermediateTarget(sigma, s, stay.particle, all.clusters, 
                                                      non.c.bar, adj, tau1, tau2, t, n, 
                                                      alpha, beta1, beta2)
  
  # if previous allocation was merge (#2), then current allocation is always merge (#2)
  if(previous.allocation == 2)
  {
    return(list("log.stay.gamma.hat" = log.stay.gamma.hat,
                "log.move.gamma.hat" = NULL,  
                "previous.allocation" = 2,
                "alternative.allocation" = NULL))
  }
  
  # otherwise split is performed (#3 or #4)
  alternative.allocation <- ifelse(particle[t-1] == 3, 4, 3)

  # determine "move" particle (i.e. choosing a different allocation decision)
  move.particle <- c(particle[1:(t-1)], alternative.allocation)
  
  # gamma hat based on move particle at time t
  log.move.gamma.hat <- LogImprovedIntermediateTarget(sigma, s, move.particle, all.clusters, 
                                                      non.c.bar, adj, tau1, tau2, t, n, 
                                                      alpha, beta1, beta2)
  
  return(list("log.stay.gamma.hat" = log.stay.gamma.hat,
              "log.move.gamma.hat" = log.move.gamma.hat,
              "previous.allocation" = previous.allocation,
              "alternative.allocation" = alternative.allocation))
}


#****************************************************
#'  Improved proposal distribution ("q hat")
#'  
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t-1 [vector]
#' @param all.clusters set of all clusters ("c") [list]
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time [scalar]
#' @param n maximum time [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter [scalar]
#' @param beta2 beta function parameter [scalar]
#' @return proposal allocation at time t
Proposal <- function(sigma, s, particle, all.clusters, non.c.bar, 
                     adj, tau1, tau2, t, n, alpha, beta1, beta2)
{
  # t=2 is a special case: 
  # allocation decisions at t=1 for all particles are #1 - so cannot use previous values
  # proposal allocation: equal chance of #2 (merge) or #4 (split)
  if(t == 2)
  {
    return(ifelse(runif(1) < 0.5, 2, 4))
  }
  
  # IF we merged at t=2 we merge throughout so that allocation = #2
  if(particle[2] == 2)
  {
    return(2)
  }
  # ELSE a split was performed at t=2: At each t there are 2 possible allocation decisions.
  # We calculate the probablity of staying in the current state (i.e. repeating the previous
  # allocation decision) by normalising over both decisions. The proposed allocation decision
  # is then made based on this probability.

  # possible allocations - and the log gamma.hats of staying in current state or moving
  possible.allocations <- PossibleAllocations(sigma, s, particle, all.clusters, non.c.bar, 
                                              adj, tau1, tau2, t, n, alpha, beta1, beta2)

  # log improved proposal probability: log of probability of staying in current state
  log.proposal.prob <- possible.allocations$log.stay.gamma.hat - 
                       logSumExp(lx = c(possible.allocations$log.stay.gamma.hat,
                                        possible.allocations$log.move.gamma.hat))
  
  # improved proposal allocation decision
  proposal.allocation <- ifelse(log(runif(1)) < log.proposal.prob, 
                                possible.allocations$previous.allocation, 
                                possible.allocations$alternative.allocation)
  
  return(proposal.allocation)
}


#****************************************************
#'  Log of improved unnormalised weights ("w hat")
#'  
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t-1 [vector]
#' @param log.previous.weight log of unnormalised weight at time t-1 [scalar]
#' @param all.clusters set of all clusters ("c") [list]
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time [scalar]
#' @param n maximum time [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter [scalar]
#' @param beta2 beta function parameter [scalar]
#' @return log unnormalised weight at time t
LogUnnormalisedWeight <- function(sigma, s, particle, log.previous.weight, all.clusters, 
                                  non.c.bar, adj, tau1, tau2, t, n, alpha, beta1, beta2)
{
  # At t=2 unnormalised weights = 2 since ratios of gamma.hats = 1 
  if(t == 2)
  {
    return(log(2)) 
  }
  
  # possible allocations - and the log gamma.hats of staying in current state or moving
  possible.allocations <- PossibleAllocations(sigma, s, particle, all.clusters, non.c.bar, 
                                              adj, tau1, tau2, t, n, alpha, beta1, beta2)
  # gamma.hat at t-1
  log.gamma.hat.previous <- LogImprovedIntermediateTarget(sigma, s, particle[1:(t-1)], 
                                                          all.clusters, non.c.bar, adj, 
                                                          tau1, tau2, t-1, n, alpha, 
                                                          beta1, beta2)
  
  # IF we merged at t=2 then we merge throughout so there is no log.move.gamma.hat
  if(particle[2] == 2)
  {
    log.weight.update <- possible.allocations$log.stay.gamma.hat - log.gamma.hat.previous
  }
  else
  {
    log.weight.update <- logSumExp(c(possible.allocations$log.stay.gamma.hat,
                                   possible.allocations$log.move.gamma.hat)) -
                         log.gamma.hat.previous
  }
  
  # update the unnormalised weight
  log.unnormalised.weight <- log.previous.weight + log.weight.update
  
  return(log.unnormalised.weight)
}

#****************************************************
#'  Resample and change history of resampled particles
#'  @param log.norm.weights log of the normalised importance weights [vector]
#'  @param particles [matrix]
#'  @param N number of particles [scalar]
#'  @param t current iteration time
#'  @return resampled particles with changed history
ResampleAndChangeParticleHistory <- function(log.norm.weights, particles, N, t)
{
  # calculate the resampling indices
  resample.indices <- StratifiedResample(log.norm.weights, N)[2:N]  # fix 1st particle  
  true.indices <- 2:N
  
  # change history of any resampled particles
  for(i in 1:(N-1))
  {
    if(true.indices[i] != resample.indices[i])
    {
      particles[true.indices[i], 1:t] <- particles[resample.indices[i], 1:t] 
    }
  }
  return(particles)
}

#****************************************************
#'  Particle Gibbs Split Merge algorithm ("Algorithm 3")
#'  
#' @param all.clusters all current clusters [list]
#' @param adj adjacency matrix for SBM [matrix]
#' @param s Anchors  [vector]
#' @param s.bar Closure of anchors [vector]
#' @param c.bar clusters that contain the anchors [list]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param N Number of particles [scalar]
#' @param resampling.threshold Resampling threshold ("beta") compared with ESS [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter 1 [scalar]
#' @param beta2 beta function parameter 2 [scalar]
#' @return Updated c.bar
ParticleGibbsSplitMerge <- function(all.clusters, adj, s, s.bar, c.bar, non.c.bar, N, 
                                    resampling.threshold, alpha, beta1, beta2)
{
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
  
  # choose a particle by multinomial sampling according to the weights
  multinomial.sample <- as.double(rmultinom(n = 1, size = 1, prob = exp(log.norm.weights)))
  updated.c.bar <- MapAllocationsToClusters(sigma, particles[which(multinomial.sample == 1), ], s)
  
  return(updated.c.bar)
}


#****************************************************
# Perform PGSMs on a selected subset of clusters
#' @param all.clusters all current clusters [list]
#' @param adj adjacency matrix for SBM [matrix]
#' @param N number of particles [scalar]
#' @param resampling.threshold resampling threshold ("beta") to be compared with ESS [scalar]
#' @param alpha tau1 parameter [scalar]
#' @param beta1 beta function parameter 1 [scalar]
#' @param beta2 beta function parameter 2 [scalar]
#' @return Updated clustering
SplitMerge <- function(all.clusters, adj, N, resampling.threshold, alpha, beta1, beta2)
{
  s <- SelectAnchors(all.clusters)  # select anchors uniformly

  # calculate c.bar and s.bar
  closure <- CalculateClosureOfAnchors(s, all.clusters)
  c.bar <- closure$c.bar
  s.bar <- closure$s.bar
  
  # calculate non.c.bar
  non.c.bar.cluster.indicators <- lapply(lapply(all.clusters, 
                                                function(x){s %in% x}), 
                                         function(y){any(y == TRUE)})
  non.c.bar <- all.clusters[which(non.c.bar.cluster.indicators == FALSE)]
  
  # update clustering
  updated.c.bar <- ParticleGibbsSplitMerge(all.clusters, adj, s, s.bar, c.bar, non.c.bar, N, 
                                           resampling.threshold, alpha, beta1, beta2)
  unchanged.clusters <- all.clusters[all.clusters %!in% c.bar]
  
  # update according to whether a split or merge was performed
  if(is.null(updated.c.bar[[2]]))
  {
    # if merge we only have 1 updated cluster
    updated.clustering <- c(list(as.integer(updated.c.bar[[1]])), 
                            unchanged.clusters)
  }
  else
  {
    updated.clustering <- c(updated.c.bar, unchanged.clusters)
  }
  
  return(updated.clustering)
}

#****************************************************
#
#************ CODE FROM RICHARD ********************* 
#
#****************************************************


#****************************************************
#' Stratified Resample
#' 
#' @param log_weights log of the normalised importance weights [vector]
#' @return resampled particles
StratifiedResample <-function(log_weights, n) # remember: n-1 since one path fixed
{
  P = length(log_weights)
  W = exp(log_weights)
  cw = cumsum(W / sum(W))
  #n = length(W) - 1
  u = runif(n,0,1)/n
  v = u + (0:(n-1))/n
  
  j = 1
  indices = matrix(0,n)
  for (i in 1:n)
  {
    while(cw[j] < v[i])
    {
      j = j+1
    }
    indices[i] = j
  }
  indices
}





logsumexp <- function(x){
  xmax = which.max(x)
  if (sum(xmax)==0)
  {
    -Inf
  }
  else
  {
    result = log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
    if (is.nan(result))
    {
      -Inf
    }
    else
    {
      result
    }
  }
}

logweightedsumexp <- function(x,w){
  logsumexp(x+w)
}

weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
}

# can take unnormalised weights
ess <- function(log_weights)
{
  #browser()
  the_ess = exp(2*logsumexp(log_weights) - logsumexp(2*log_weights))
  
  if (is.nan(the_ess))
  {
    #browser()
    0#length(log_weights)
    #normlogweights = log_0weights-logsumexp(log_weights)
    #1.0/sum(exp(2.0*log_weights))
  }
  else
  {
    the_ess
  }
}

