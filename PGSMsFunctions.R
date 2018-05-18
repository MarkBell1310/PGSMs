
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
#'  Pre-compute edge counts BETWEEN c.bar and non.c.bar clusters
#'  
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param directed whether network is directed [boolean]
#' 
#' @return list of edge "counts" [list of vectors] 
PreComputeEdgeCountsBetweenCbarAndNonCbar <- function(adj, non.c.bar, sigma, directed)
{
  # returns list of length(non.c.bar)
  # each vector in list is of length(sigma):
  # -the first n = length(c.bar[[1]]) elements are counts with c.bar cluster 1
  # -next m = length(c.bar[[2]]) elements are counts with c.bar cluster 2 (if it exists)
  
  # count edges between each node in sigma (ordered nodes of c.bar) & each cluster of non.c.bar 
  if(directed == FALSE)
  {
    return(lapply(non.c.bar, function(x)
    {
      sapply(sigma, function(z)
      {
        sum(adj[z, x])
      })
    }))
  }
  
  if(directed == TRUE)
  {
    return(lapply(non.c.bar, function(x)
    {
      sapply(sigma, function(z)
      {
        sum(adj[z, x]) + sum(adj[x, z])
      })
    }))
  }
}
#PreComputeEdgeCountsBetweenCbarAndNonCbar(adj, non.c.bar, sigma, directed)
#pre.computed.edge.counts <- PreComputeEdgeCountsBetweenCbarAndNonCbar(adj, non.c.bar, sigma, directed)

#****************************************************
#'  Count new edges WITHIN the c.bar clusters 
#'  (between the new node and remaining nodes)
#'  
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' 
#' @return new edge counts WITHIN the c.bar clusters [vector] 
#'         [vector length = 2]
CountNewEdgesWithinCbarClusters <- function(c.bar.current, adj, t, sigma, particle, directed)
{
  ## Count edges between new node and all remaining nodes in the same cluster
  ## (that new node has been added to)
  
  ## Output given in vector of form below - one of the two values must be zero:
  ## c(within.c.bar1, within.c.bar2)
  
  # for undirected networks only need to sum adjacency matrix one way
  if(directed == FALSE)
  {
    # IF only 1 cluster
    # OR if 2 clusters AND new node added to cluster 1
    if(length(c.bar.current) == 1 || particle[t] == 3)
    {
      return(c(sum(adj[sigma[t], c.bar.current[[1]]]), 0))
      #max.counts <- c(length(c.bar.current[[1]]) - 1, 0)
    }
    
    # ELSE if 2 clusters AND new node added to cluster 2 
    if(particle[t] == 4)
    {
      return(c(0, sum(adj[sigma[t], c.bar.current[[2]]])))
      #max.counts <- c(0, length(c.bar.current[[2]]) - 1)
    }
  }

  # for directed networks need to sum adjacency matrix both ways
  if(directed == TRUE)
  {
    # IF only 1 cluster
    # OR if 2 clusters AND new node added to cluster 1
    if(length(c.bar.current) == 1 || particle[t] == 3)
    {
      return(c(sum(adj[sigma[t], c.bar.current[[1]]]) + 
                 sum(adj[c.bar.current[[1]], sigma[t]]), 0))
      #max.counts <- c(2 * (length(c.bar.current[[1]]) - 1), 0) 
    }
    
    # ELSE if 2 clusters AND new node added to cluster 2 
    if(particle[t] == 4)
    {
      return(c(0, sum(adj[sigma[t], c.bar.current[[2]]]) + 
                    sum(adj[c.bar.current[[2]], sigma[t]])))
      #max.counts <- c(0, 2 * (length(c.bar.current[[2]]) - 1)) 
    }
  }
}
#CountNewEdgesWithinCbarClusters(c.bar.current, adj, t, sigma, particle, directed = FALSE)


#****************************************************
#'  Count new edges BETWEEN the c.bar clusters
#'  
#' @param c.bar.current The 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' 
#' @return new edge counts BETWEEN the c.bar clusters [scalar] 
CountNewEdgesBetweenCbarClusters <- function(c.bar.current, adj, t, sigma, particle, directed)
{
  # count edges between new node and all elements in the OTHER cluster
  
  # for undirected networks only need to sum adjacency matrix one way
  if(directed == FALSE)
  {
    # IF new node added to cluster 2
    if(particle[t] == 4)
    {
      return(sum(adj[sigma[t], c.bar.current[[1]]]))
      #max.counts <- length(c.bar.current[[1]])
    }
    
    # IF new node added to cluster 1
    if(particle[t] == 3)
    {
      return(sum(adj[sigma[t], c.bar.current[[2]]]))
      #max.counts <- length(c.bar.current[[2]])
    }
  }
  
  # for directed networks need to sum adjacency matrix both ways
  if(directed == TRUE)
  {
    # IF new node added to cluster 2
    if(particle[t] == 4)
    {
      return(sum(adj[sigma[t], c.bar.current[[1]]]) +
             sum(adj[c.bar.current[[1]], sigma[t]]))
      #max.counts <- length(c.bar.current[[1]]) * 2
    }
    
    # IF new node added to cluster 1
    if(particle[t] == 3)
    {
      return(sum(adj[sigma[t], c.bar.current[[2]]]) +
             sum(adj[c.bar.current[[2]], sigma[t]]))
      #max.counts <- length(c.bar.current[[2]]) * 2
    }
  }
}
#CountNewEdgesBetweenCbarClusters(c.bar.current, adj, t, sigma, particle, directed = FALSE)


#****************************************************
#'  Count new edges BETWEEN c.bar and non.c.bar clusters
#'  
#' @param pre.computed.edge.counts The pre-computed edge counts [list of vectors]
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' 
#' @return new edge counts BETWEEN c.bar and non.c.bar clusters [vector]
#'         [vector length = 2 * length(non.c.bar)]
CountNewEdgesBetweenCbarAndNonCbar <- function(pre.computed.edge.counts, c.bar.current, 
                                               t, particle)
{
  ## Determine which c.bar cluster the new node was added to, then count edges but use
  ## correct form for running totals vectors - which means leaving certain entries zero.
  
  # Node added to cluster 1:
  # IF only 1 cluster OR if 2 clusters AND new node added to cluster 1
  if(length(c.bar.current) == 1 || particle[t] == 3)
  {
    return(c(sapply(pre.computed.edge.counts, function(x){x[t]}),
                rep(0, length(pre.computed.edge.counts))))
  }
  
  # Node added to cluster 2:
  # ELSE if 2 clusters AND new node added to cluster 2 
  if(particle[t] == 4)
  {
    return(c(rep(0, length(pre.computed.edge.counts)),
                sapply(pre.computed.edge.counts, function(x){x[t]})))
  }
}
#CountNewEdgesBetweenCbarAndNonCbar(pre.computed.edge.counts, c.bar.current, t, particle)

#****************************************************
#'  Count new edges for all 3 types
#'  
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' 
#' @return new edge counts for all 3 types; in the correct form for RunningTotalEdgeCount()
#'         function  [vector]
NewEdgeCounts <- function(c.bar.current, adj, non.c.bar, sigma, particle, t, directed)  
{
  # This function needs to consider the following 3 types of edge counts:
  # (1) nodes WITHIN cluster(s) in c.bar.current
  # (2) nodes BETWEEN the 2 clusters in c.bar.current (only if we have 2 clusters) 
  # (3) nodes BETWEEN c.bar cluster(s) and non.c.bar cluster(s)
  
  # Output vector assumes the following (2 clusters in c.bar) form:
  # c(WithinCbar1, WithinCbar2, BetweenCbar, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...
  #   BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...)
  
  # TYPE (1): nodes within clusters in c.bar
  within.c.bar.current <- CountNewEdgesWithinCbarClusters(c.bar.current, adj, t, sigma, 
                                                          particle, directed)
  
  # TYPE (2): nodes between clusters in c.bar (c.bar must have 2 clusters)
  if(length(c.bar.current) == 2)
  {
    between.c.bar.current <- CountNewEdgesBetweenCbarClusters(c.bar.current, adj, t, sigma, 
                                                              particle, directed)
  }
  else
  {
    between.c.bar.current <- 0
  }
  
  # TYPE (3): nodes between c.bar and non c.bar clusters
  if(length(non.c.bar) != 0)
  {
    between.c.bar.non.c.bar <- CountNewEdgesBetweenCbarAndNonCbar(pre.computed.edge.counts, 
                                                                  c.bar.current, t, particle)
  }
  else
  {
    between.c.bar.non.c.bar <- rep(0, 2 * length(non.c.bar))
  }
  
  return(c(within.c.bar.current, between.c.bar.current, between.c.bar.non.c.bar))
}
#NewEdgeCounts(c.bar.current, adj, non.c.bar, sigma, particle, t, directed = FALSE)  


#****************************************************
#'  Calculate maximum edge counts 
#'  
#' @param previous.total Previous running totals of edge counts [vector]
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' 
#' @return maximum edge counts [vector]
MaximumEdgeCounts <- function(previous.total, c.bar.current, adj, non.c.bar, sigma, 
                              particle, t, directed)
{
  # TO DO: inlcude if statements directed network
  # TO DO: create vector in same form as for edge counts
  
  if(directed == TRUE)
  {
    
  }
  
  if(directed == FALSE)
  {
    
  }
  
  #### Code for calculating max edges for type (3) counts 
  # # max edge counts BETWEEN all c.bar.current and non.c.bar clusters
  # max.counts.between.c.bar.non.c.bar.clusters <- as.double(sapply(non.c.bar, function(x)
  # {
  #   if(length(c.bar.current) == 1) # if c.bar.current only contains 1 cluster
  #   {
  #     return(length(c.bar.current[[1]]) * length(x))
  #   }
  #   if(length(c.bar.current) == 2)
  #   {
  #     return(c(length(c.bar.current[[1]]) * length(x), length(c.bar.current[[2]]) * length(x)))
  #   }
  # }))
  
  # max.counts <- as.double(sapply(non.c.bar, function(x)
  # {
  #   if(length(c.bar.current) == 1) 
  #   {
  #     return(length(c.bar.current[[1]]) * length(x))
  #   }
  #   if(length(c.bar.current) == 2)
  #   {
  #     return(c(length(c.bar.current[[1]]) * length(x), length(c.bar.current[[2]]) * length(x)))
  #   }
  # }))
  
  
  #*************************
  # TO DO: consider - can we still use the same general (edge count) form if only have 1 cluster??
  # How to remove "redundant" running totals for likelihood? Since some non-redundant ones will be zero.
  # How to calculate max.counts? Don't calculate at each time update - do this in new EdgeCounts() function
  #*************************
  
}

#****************************************************
#'  Calculate total edge counts 
#'  
#' @param previous.total Previous running totals of edge counts [vector]
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' 
#' @return total edge counts [vector]
TotalEdgeCounts <- function(previous.total, c.bar.current, adj, non.c.bar, sigma, 
                            particle, t, directed)
{
  # TO DO: pass the previous.total down from higher level functions
  
  ### 2 "FORMS" of vector for storing edge counts. So far we have assumed FORM 1.
  ### Now need to convert to appropriate form - depending on no. c.bar clusters 
  # "FORM 1": 2 clusters in c.bar (2m + 3 elements), where m = no. non.c.bar clusters
  #   c(WithinCbar1, WithinCbar2, BetweenCbar, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...
  #     BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...)
  # "FORM 2": 1 cluster in c.bar (m + 1 elements)
  #   c(WithinCbar1, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...)
  
  
  # IF split particle then leave new.total in its current form below ("FORM 1"):
  #   c(WithinCbar1, WithinCbar2, BetweenCbar, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ... ,
  #     BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...)
  new.total <- previous.total + NewEdgeCounts(c.bar.current, adj, non.c.bar, sigma, 
                                              particle, t, directed) 

  # IF merge particle then remove redundant elements of new.total ("FORM 2"):
  #   c(WithinCbar1, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...)
  if(particle[2] == 2)
  {
    m <- length(non.c.bar)
    new.total <- new.total[-c(2, 3, (3 + m + 1):length(new.total))]
  }
  
  return(new.total)
}


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

