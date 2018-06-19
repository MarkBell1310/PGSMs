
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

  # 1st anchor point (from sigma) goes in cluster 1
  cluster1[1] <- sigma[1] 
  cluster2[1] <- sigma[2]
  
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
#'  Count new edges WITHIN the c.bar clusters 
#'  (between the new node and remaining nodes)
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
#' @return new edge counts WITHIN the c.bar clusters [vector: length = 2] 
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
      return(list("edge.counts" = c(sum(adj[sigma[t], c.bar.current[[1]]]), 0),
                  "max.counts" = c(length(c.bar.current[[1]]) - 1, 0)))
    }
    
    # ELSE if 2 clusters AND new node added to cluster 2 
    if(particle[t] == 4)
    {
      return(list("edge.counts" = c(0, sum(adj[sigma[t], c.bar.current[[2]]])),
                  "max.counts" = c(0, length(c.bar.current[[2]]) - 1)))
    }
  }

  # for directed networks need to sum adjacency matrix both ways
  if(directed == TRUE)
  {
    # IF only 1 cluster
    # OR if 2 clusters AND new node added to cluster 1
    if(length(c.bar.current) == 1 || particle[t] == 3)
    {
      return(list("edge.counts" = c(sum(adj[sigma[t], c.bar.current[[1]]]) + 
                                      sum(adj[c.bar.current[[1]], sigma[t]]), 0),
                  "max.counts" = c(2 * (length(c.bar.current[[1]]) - 1), 0)))
    }
    
    # ELSE if 2 clusters AND new node added to cluster 2 
    if(particle[t] == 4)
    {
      return(list("edge.counts" = c(0, sum(adj[sigma[t], c.bar.current[[2]]]) + 
                                      sum(adj[c.bar.current[[2]], sigma[t]])),
                  "max.counts" = c(0, 2 * (length(c.bar.current[[2]]) - 1)))) 
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
      return(list("edge.counts" = sum(adj[sigma[t], c.bar.current[[1]]]),
                  "max.counts" = length(c.bar.current[[1]])))
    }
    
    # IF new node added to cluster 1
    if(particle[t] == 3)
    {
      return(list("edge.counts" = sum(adj[sigma[t], c.bar.current[[2]]]),
                  "max.counts" = length(c.bar.current[[2]])))
    }
  }
  
  # for directed networks need to sum adjacency matrix both ways
  if(directed == TRUE)
  {
    # IF new node added to cluster 2
    if(particle[t] == 4)
    {
      return(list("edge.counts" = sum(adj[sigma[t], c.bar.current[[1]]]) +
                                  sum(adj[c.bar.current[[1]], sigma[t]]),
                  "max.counts" = length(c.bar.current[[1]]) * 2))
    }
    
    # IF new node added to cluster 1
    if(particle[t] == 3)
    {
      return(list("edge.counts" = sum(adj[sigma[t], c.bar.current[[2]]]) +
                                  sum(adj[c.bar.current[[2]], sigma[t]]),
                  "max.counts" = length(c.bar.current[[2]]) * 2))
    }
  }
}
#CountNewEdgesBetweenCbarClusters(c.bar.current, adj, t, sigma, particle, directed = FALSE)


#****************************************************
#'  Count new edges BETWEEN c.bar and non.c.bar clusters
#' @param global.counts.between.nodes.clusters Edge counts between every node & cluster [matrix]
#' @param global.non.c.bar.indices Indices of non.c.bar clusters [vector]
#' @param c.bar.current The 1 or 2 clusters that contain the anchors [list]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param sigma Ordered elements of c.bar [vector]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param directed whether network is directed [boolean]
#' @return new edge counts BETWEEN c.bar and non.c.bar clusters [vector]
#'         [vector length = 2 * length(non.c.bar)]
CountNewEdgesBetweenCbarAndNonCbar <- function(global.counts.between.nodes.clusters, 
                                               global.non.c.bar.indices, c.bar.current, 
                                               t, particle, sigma, non.c.bar, directed)
{
  # Count edges between new node & each cluster of non.c.bar
  new.edge.counts <- 
    global.counts.between.nodes.clusters[sigma[t], global.non.c.bar.indices] 
  m <- length(new.edge.counts)
  
  # if network undirected, max counts is correct - otherwise double the values
  max.counts <- sapply(non.c.bar, length)
  
  if(directed == TRUE)
  {
    max.counts <- 2 * max.counts
                         
  }
  
  # Node added to cluster 1:
  # IF only 1 cluster OR if 2 clusters AND new node added to cluster 1
  if(length(c.bar.current) == 1 || particle[t] == 3)
  {
    return(list("edge.counts" = c(new.edge.counts, rep(0, m)),
                "max.counts" = c(max.counts, rep(0, m))))
  }
  
  # Node added to cluster 2:
  # ELSE if 2 clusters AND new node added to cluster 2 
  if(particle[t] == 4)
  {
    return(list("edge.counts" = c(rep(0, m), new.edge.counts),
                "max.counts" = c(rep(0, m), max.counts)))
  }
}
#CountNewEdgesBetweenCbarAndNonCbar(global.counts.between.nodes.clusters, global.non.c.bar.indices, c.bar.current, t, particle, sigma, non.c.bar, directed)

#****************************************************
#'  Count new edges for all 3 types
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed whether network is directed [boolean]
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
  within.c.bar <- CountNewEdgesWithinCbarClusters(c.bar.current, adj, t, sigma, 
                                                  particle, directed)
  within.c.bar.edge.counts <- within.c.bar$edge.counts
  within.c.bar.max.counts <- within.c.bar$max.counts
  
  # TYPE (2): nodes between clusters in c.bar (c.bar must have 2 clusters)
  if(length(c.bar.current) == 2)
  {
    between.c.bar <- CountNewEdgesBetweenCbarClusters(c.bar.current, adj, t, sigma, 
                                                      particle, directed)
    between.c.bar.edge.counts <- between.c.bar$edge.counts
    between.c.bar.max.counts <- between.c.bar$max.counts
  }
  else
  {
    between.c.bar.edge.counts <- between.c.bar.max.counts <- 0
  }
  
  # TYPE (3): nodes between c.bar and non c.bar clusters
  if(length(non.c.bar) >= 1)
  {
    between.c.bar.non.c.bar <- 
      CountNewEdgesBetweenCbarAndNonCbar(global.counts.between.nodes.clusters, 
                                         global.non.c.bar.indices, c.bar.current, 
                                         t, particle, sigma, non.c.bar, directed)
    between.c.bar.non.c.bar.edge.counts <- between.c.bar.non.c.bar$edge.counts
    between.c.bar.non.c.bar.max.counts <- between.c.bar.non.c.bar$max.counts
  }
  else
  {
    # if no non.c.bar clusters, set to NULL
    between.c.bar.non.c.bar.edge.counts <- between.c.bar.non.c.bar.max.counts <- NULL
  }
  
  return(list("edge.counts" = c(within.c.bar.edge.counts, between.c.bar.edge.counts, 
                                between.c.bar.non.c.bar.edge.counts),
              "max.counts" = c(within.c.bar.max.counts, between.c.bar.max.counts, 
                               between.c.bar.non.c.bar.max.counts)))
}
#NewEdgeCounts(c.bar.current, adj, non.c.bar, sigma, particle, t, directed = FALSE)  


#****************************************************
#'  Calculate running total of edge counts 
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return total edge counts [vector]
RunningTotalsEdgeCounts <- function(c.bar.current, adj, non.c.bar, sigma, particle, t, directed,
                                  particle.index)
{
  # new edge counts
  new.counts <- NewEdgeCounts(c.bar.current, adj, non.c.bar, sigma, particle, t, directed) 
  new.edge.counts <- new.counts$edge.counts
  new.max.counts <- new.counts$max.counts
  
  # update running totals - but don't update global variables until particle chosen
  edge.count.running.total <- global.running.total.edge.counts[[particle.index]] + new.edge.counts
  max.count.running.total <- global.running.total.max.counts[[particle.index]] + new.max.counts

  return(list("edge.count.running.total" = edge.count.running.total,
              "max.count.running.total" = max.count.running.total))
}
#RunningTotalsEdgeCounts(c.bar.current, adj, non.c.bar, sigma, particle, t, directed,particle.index)

#****************************************************
#'  Counts for likelihood - Remove redundant elements of count vectors
#' @param c.bar.current clusters that contain the anchors, filled in up to time t [list]
#' @param adj adjacency matrix of the SBM [matrix]
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param t current time [scalar]
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return total edge counts and maximum edge counts [list of vectors]
CountsForLikelihood <- function(c.bar.current, adj, non.c.bar, sigma, particle, 
                                t, directed, particle.index)
{
  ### Note: 2 "FORMS" of vector for storing edge counts. So far we have assumed FORM 1.
  ### Now need to convert to appropriate form - depending on no. c.bar clusters 
  # "FORM 1": 2 clusters in c.bar (2m + 3 elements), where m = no. non.c.bar clusters
  #   c(WithinCbar1, WithinCbar2, BetweenCbar, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...
  #     BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...)
  # "FORM 2": 1 cluster in c.bar (m + 1 elements)
  #   c(WithinCbar1, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...)
  
  ### Note: if no non.c.bar clusters, then terms (BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...
  #         BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...) are NULL. 
  
  ##### Split particle
  # IF split particle then leave edge counts in their current form below ("FORM 1"):
  #   c(WithinCbar1, WithinCbar2, BetweenCbar, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ... ,
  #     BetweenCbar2NonCbar1, BetweenCbar2NonCbar2, ...)
  running.totals <- RunningTotalsEdgeCounts(c.bar.current, adj, non.c.bar, sigma, 
                                            particle, t, directed, particle.index)
  edge.count.running.total <- running.totals$edge.count.running.total
  max.count.running.total <- running.totals$max.count.running.total
  
  likelihood.edge.count.total <- edge.count.running.total
  likelihood.max.count.total <- max.count.running.total
  
  ##### Merge particle
  # IF merge particle then remove redundant elements of vectors for likelihood calc ("FORM 2"):
  #   c(WithinCbar1, BetweenCbar1NonCbar1, BetweenCbar1NonCbar2, ...)
  if(particle[2] == 2)
  {
    m <- length(non.c.bar)
    
    # IF there is at least one non.c.bar cluster
    if(m > 0)
    {
      likelihood.edge.count.total <- 
        edge.count.running.total[-c(2, 3, (3 + m + 1):length(edge.count.running.total))]
      likelihood.max.count.total <- 
        max.count.running.total[-c(2, 3, (3 + m + 1):length(max.count.running.total))] 
    }
    
    # IF there are no non.c.bar clusters
    if(m == 0)
    {
      likelihood.edge.count.total <- edge.count.running.total[-c(2, 3)]
      likelihood.max.count.total <- max.count.running.total[-c(2, 3)]
    }  
  }
  return(list("edge.count.running.total" = edge.count.running.total,
              "max.count.running.total" = max.count.running.total,
              "likelihood.edge.count.total" = likelihood.edge.count.total,
              "likelihood.max.count.total" = likelihood.max.count.total))
}
#CountsForLikelihood(c.bar.current, adj, non.c.bar, sigma, particle, t, directed, particle.index)

#****************************************************
#'  Log of intermediate target distribution at time t (Log "Gamma")
#'
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param all.clusters set of all clusters ("c") [list]
#' @param non.c.bar clusters not containing anchors [list]
#' @param adj adjacency matrix of SBM [matrix]
#' @param tau1 factorisation of prior [function]
#' @param tau2 factorisation of prior [function]
#' @param t current time
#' @param alpha tau1 parameter
#' @param beta1 beta function parameter
#' @param beta2 beta function parameter
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return log of intermeduate target distribution
LogIntermediateTarget <- function(sigma, s, particle, all.clusters, non.c.bar, adj, 
                                  tau1, tau2, t, alpha, beta1, beta2, directed, 
                                  particle.index)
{
  # calculate "c.bar.current": c.bar at time t
  c.bar.current <- MapAllocationsToClusters(sigma[1:t], particle, s)
  if(is.null(c.bar.current[[2]]))
  {
    c.bar.current <- list(c.bar.current[[1]])
  }
  
  # count relevant edges within and between clusters
  all.counts <- CountsForLikelihood(c.bar.current, adj, non.c.bar, sigma, particle, 
                                    t, directed, particle.index)
  edge.counts <- all.counts$likelihood.edge.count.total
  max.counts <- all.counts$likelihood.max.count.total
  n <- length(edge.counts)
  
  # product of all between-cluster log likelihoods 
  log.likelihoods <- rep(0, n)
  for (i in 1:n)
  {
    log.likelihoods[i] <-  log(beta(beta1 + edge.counts[i], 
                                    max.counts[i] - edge.counts[i] + beta2)) - 
                            log(beta(beta1, beta2))
  }
  
  sum.log.likelihoods <- sum(log.likelihoods)
  
  # prior (3rd line is 0 if we split)
  log.prior <- log(tau1(alpha, length(c.bar.current) + length(all.clusters) - length(c.bar.current))) +
               log(tau2(length(c.bar.current[[1]]))) + 
               ifelse(particle[t] != 2, log(tau2(length(c.bar.current[[2]]))), 0) 

  # intermediate target
  log.int.target <- log.prior + sum.log.likelihoods
  
  return(list("log.int.target" = log.int.target,
              "edge.count.running.total" = all.counts$edge.count.running.total,
              "max.count.running.total" = all.counts$max.count.running.total))
}   
#LogIntermediateTarget(sigma, s, particle, all.clusters, c.bar.current, non.c.bar, adj, tau1, tau2, t, alpha, beta1, beta2, directed, particle.index)

#****************************************************
#' Log of Improved intermediate target distribution (Log "Gamma hat")
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
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return log of improved intermeduate target distribution
LogImprovedIntermediateTarget <- function(sigma, s, particle, all.clusters, non.c.bar, 
                                          adj, tau1, tau2, t, n, alpha, beta1, beta2, 
                                          directed, particle.index)
{
  # log intermediate target at current time t
  log.intermediate.target <- LogIntermediateTarget(sigma, s, particle, all.clusters, 
                                                   non.c.bar, adj, tau1, tau2, t, alpha, 
                                                   beta1, beta2, directed, particle.index)
  log.gamma_t <- log.intermediate.target$log.int.target
  
  # log improved intermediate target ("gamma hat"): also uses gamma at t=2
  log.improved.int.target <- (zeta_t(t, n) - 1) * global.log.gamma_2[particle.index] + 
                              log.gamma_t
  
  return(list("log.imp.int.target" = log.improved.int.target,
              "edge.count.running.total" = log.intermediate.target$edge.count.running.total,
              "max.count.running.total" = log.intermediate.target$max.count.running.total))
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
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return the possible allocations and resulting log "gamma hats" at time t
PossibleAllocations <- function(sigma, s, particle, all.clusters, non.c.bar, adj, tau1, tau2, 
                                t, n, alpha, beta1, beta2, directed, particle.index)
{ 
  previous.allocation <- particle[t-1]
  
  # particle and log gamma.hat resulting from staying in current state 
  # (i.e. repeating previous allocation decision)
  stay.particle <- c(particle[1:(t-1)], previous.allocation)
  log.stay.improved.int.target <- 
    LogImprovedIntermediateTarget(sigma, s, stay.particle, all.clusters, non.c.bar, adj, 
                                  tau1, tau2, t, n, alpha, beta1, beta2, directed, particle.index)
  log.stay.gamma.hat <- log.stay.improved.int.target$log.imp.int.target
  stay.edge.counts <- log.stay.improved.int.target$edge.count.running.total
  stay.max.counts <- log.stay.improved.int.target$max.count.running.total
  
  # if previous allocation was merge (#2), then current allocation is always merge (#2)
  if(previous.allocation == 2)
  {
    return(list("log.stay.gamma.hat" = log.stay.gamma.hat,
                "log.move.gamma.hat" = NULL,  
                "previous.allocation" = 2,
                "alternative.allocation" = NULL,
                "stay.edge.counts" = stay.edge.counts,
                "stay.max.counts" = stay.max.counts))
  }
  
  # otherwise split is performed (#3 or #4)
  alternative.allocation <- ifelse(particle[t-1] == 3, 4, 3)

  # determine "move" particle (i.e. choosing a different allocation decision)
  move.particle <- c(particle[1:(t-1)], alternative.allocation)
  
  # gamma hat based on move particle at time t
  log.move.improved.int.target <- 
    LogImprovedIntermediateTarget(sigma, s, move.particle, all.clusters, non.c.bar, adj, 
                                  tau1, tau2, t, n, alpha, beta1, beta2, directed, particle.index)
  log.move.gamma.hat <- log.move.improved.int.target$log.imp.int.target
  move.edge.counts <- log.move.improved.int.target$edge.count.running.total
  move.max.counts <- log.move.improved.int.target$max.count.running.total
  
  return(list("log.stay.gamma.hat" = log.stay.gamma.hat,
              "log.move.gamma.hat" = log.move.gamma.hat,
              "previous.allocation" = previous.allocation,
              "alternative.allocation" = alternative.allocation,
              "stay.edge.counts" = stay.edge.counts,
              "move.edge.counts" = move.edge.counts,
              "stay.max.counts" = stay.max.counts,
              "move.max.counts" = move.max.counts))
}
#PossibleAllocations(sigma, s, particle, all.clusters, non.c.bar, adj, tau1, tau2, t, n, alpha, beta1, beta2, directed, particle.index)


#****************************************************
#'  Improved proposal distribution ("q hat")
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
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @return proposal allocation at time t
Proposal <- function(sigma, s, particle, all.clusters, non.c.bar, adj, tau1, tau2, 
                     t, n, alpha, beta1, beta2, directed, particle.index)
{
  ##--------------------------------------
  ## CASE 1: t = 2 
  # log gamma hats = 0 but need to store edge counts for running totals and global.log.gamma_2
  # allocation decisions at t=1 for all particles are #1 - so cannot use previous values
  # proposal allocation (for p>1): equal chance of decision #2 (merge) or #4 (split)
  if(t == 2)
  { 
    # determine chosen allocation decision and resulting particle
    allocation <- ifelse(runif(1) < 0.5, 2, 4)
    if(particle.index > 1)
    {
      particle <- c(1, allocation )
    }
    
    # calculate edge count running totals and global_log.gamma_2
    # (global.log.gamma_2 not used until t = 3)
    log.intermediate.target <-  
      LogIntermediateTarget(sigma, s, particle, all.clusters, non.c.bar, adj, tau1, tau2, 
                            t, alpha, beta1, beta2, directed, particle.index)
    
    # update edge count running totals
    global.running.total.edge.counts[[particle.index]] <<- 
      log.intermediate.target$edge.count.running.total
    global.running.total.max.counts[[particle.index]] <<- 
      log.intermediate.target$max.count.running.total
    
    # store global.log.gamma_2
    global.log.gamma_2[particle.index] <<- log.intermediate.target$log.int.target
    
    return(list("allocation" = allocation,
                "log.stay.gamma.hat" = 0,
                "log.move.gamma.hat" = 0,
                "log.gamma.hat" = 0))
  }
  
  ##--------------------------------------
  ## CASE 2: t > 2, any merge particles
  ## -only needs the "stay" scenario
  
  # possible allocations - and the log gamma.hats of staying in current state or moving
  possible.allocations <- PossibleAllocations(sigma, s, particle, all.clusters, non.c.bar, 
                                              adj, tau1, tau2, t, n, alpha, beta1, beta2,
                                              directed, particle.index)
  
  # IF we merged at t=2 we merge throughout so that allocation = #2
  if(particle[2] == 2)
  {
    global.running.total.edge.counts[[particle.index]] <<- possible.allocations$stay.edge.counts
    global.running.total.max.counts[[particle.index]] <<- possible.allocations$stay.max.counts
    
    return(list("allocation" = 2,
                "log.stay.gamma.hat" = possible.allocations$log.stay.gamma.hat,
                "log.move.gamma.hat" = NULL, 
                "log.gamma.hat" = possible.allocations$log.stay.gamma.hat))
  }
  
  # ELSE a split was performed at t=2: At each t there are 2 possible allocation decisions.
  # We calculate the probablity of staying in the current state (i.e. repeating the previous
  # allocation decision) by normalising over both decisions. The proposed allocation decision
  # is then made based on this probability.

  # log improved proposal probability: log of probability of staying in current state
  log.proposal.prob <- possible.allocations$log.stay.gamma.hat - 
                       logSumExp(lx = c(possible.allocations$log.stay.gamma.hat,
                                        possible.allocations$log.move.gamma.hat))
  
  ##--------------------------------------
  ## CASE 3: t > 2, particle.index = 1 (where 1st particle is a split)
  #  -Use particle from conditional path to calculate log gamma hat
  #  -also store the corresponding edge count running totals
  #  -still need to output the log.stay.gamma.hat and log.move.gamma.hat - for the weights  
  if(particle.index == 1 && t > 2)
  {
    # calculate edge count running total and log gamma hat
    log.imp.int.target <-  
      LogImprovedIntermediateTarget(sigma, s, particle, all.clusters, non.c.bar, 
                                    adj, tau1, tau2, t, n, alpha, beta1, beta2, 
                                    directed, particle.index)
    
    # update edge count running totals
    global.running.total.edge.counts[[particle.index]] <<- 
      log.imp.int.target$edge.count.running.total
    global.running.total.max.counts[[particle.index]] <<- 
      log.imp.int.target$max.count.running.total
    
    return(list("allocation" = particle[t],
                "log.stay.gamma.hat" = possible.allocations$log.stay.gamma.hat,
                "log.move.gamma.hat" = possible.allocations$log.move.gamma.hat,
                "log.gamma.hat" = log.imp.int.target$log.imp.int.target))
  }
  
  ##--------------------------------------
  ## CASE 4: t > 2, any split particles with particle.index > 1 
  #  -Select proposal allocation decision and corresponding log.gamma.hat 
  #  -Edge counts added to running totals
  
  # IF we "stay" (repeat) previous decision
  if(log(runif(1)) < log.proposal.prob)
  {
    proposal.allocation <- possible.allocations$previous.allocation
    log.gamma.hat <- possible.allocations$log.stay.gamma.hat
    global.running.total.edge.counts[[particle.index]] <<- possible.allocations$stay.edge.counts
    global.running.total.max.counts[[particle.index]] <<- possible.allocations$stay.max.counts
  }
  # ELSE we "move" (choose different) decision
  else
  {
    proposal.allocation <- possible.allocations$alternative.allocation
    log.gamma.hat <- possible.allocations$log.move.gamma.hat
    global.running.total.edge.counts[[particle.index]] <<- possible.allocations$move.edge.counts
    global.running.total.max.counts[[particle.index]] <<- possible.allocations$move.max.counts
  }

  
  return(list("allocation" = proposal.allocation,
              "log.stay.gamma.hat" = possible.allocations$log.stay.gamma.hat,
              "log.move.gamma.hat" = possible.allocations$log.move.gamma.hat,
              "log.gamma.hat" = log.gamma.hat))
}


#****************************************************
#'  Log of improved unnormalised weights ("w hat")
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t-1 [vector]
#' @param log.previous.unnormalised.weight log of unnormalised weight at time t-1 [scalar]
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
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @param log.gamma.hat.previous log intermediate target at time t-1 [scalar]
#' @return log unnormalised weight at time t
LogUnnormalisedWeight <- function(sigma, s, particle, log.previous.unnormalised.weight, 
                                  all.clusters, non.c.bar, adj, tau1, tau2, t, n, alpha, 
                                  beta1, beta2, directed, particle.index, log.gamma.hat.previous)
{
  # Proposal 
  proposal <- Proposal(sigma, s, particle, all.clusters, non.c.bar, adj, tau1, tau2, t, n, 
                       alpha, beta1, beta2, directed, particle.index)

  # At t=2 unnormalised weights = 2 since ratios of gamma.hats = 1 
  if(t == 2)
  {
    return(list("log.unnormalised.weight" = log(2),
                "proposal.allocation" = proposal$allocation,
                "log.gamma.hat" = 0))
  }

  log.stay.gamma.hat <- proposal$log.stay.gamma.hat
  log.move.gamma.hat <- proposal$log.move.gamma.hat
  
  # IF we merged at t=2 then we merge throughout so there is no log.move.gamma.hat
  if(particle[2] == 2)
  {
    log.weight.update <- log.stay.gamma.hat - log.gamma.hat.previous
  }
  else
  {
    log.weight.update <- logSumExp(c(log.stay.gamma.hat, log.move.gamma.hat)) -
                          log.gamma.hat.previous
  }
  
  # update the unnormalised weight
  log.unnormalised.weight <- log.previous.unnormalised.weight + log.weight.update
  
  return(list("log.unnormalised.weight" = log.unnormalised.weight,
              "proposal.allocation" = proposal$allocation,
              "log.gamma.hat" = proposal$log.gamma.hat))
}


#****************************************************
#'  Create global matrix of edge counts between every node and every cluster
#' @param all.clusters all current clusters [list of vectors]
#' @param nun.nodes number of nodes [scalar]
#' @param directed whether network is directed or not [boolean]
#' @return edge counts [matrix: dim = no. nodes x no. clusters]
CreateGlobalMatrixEdgeCountsBetweenAllNodesAndClusters <- function(all.clusters, num.nodes, 
                                                                   directed)
{
  K <- length(all.clusters)
  edge.counts.matrix <- Matrix(rep(0, num.nodes*K), c(num.nodes, K))
  
  # loop through all nodes and all clusters
  for(node in 1:num.nodes)
  {
    for(cluster in 1:K)
    {
      # if network undirected
      if(directed == FALSE)
      {
        edge.counts.matrix[node, cluster] <- sum(adj[node, all.clusters[[cluster]]])
      }
      
      # if network directed
      if(directed == TRUE)
      {
        edge.counts.matrix[node, cluster] <- sum(adj[node, all.clusters[[cluster]]]) + 
                                              sum(adj[all.clusters[[cluster]], node])
      }
    }
  }
  return(edge.counts.matrix)
}
#CreateGlobalMatrixEdgeCountsBetweenAllNodesAndClusters(all.clusters, num.nodes, directed = FALSE)

#****************************************************
#'  Create global matrix of edge counts between every node and every cluster
#' @param all.clusters all current clusters [list of vectors]
#' @param nun.nodes number of nodes [scalar]
#' @return edge counts [matrix: dim = no. clusters x no. clusters]
CreateGlobalMatrixEdgeCountsBetweenClusters <- function(all.clusters, directed)
{
  K <- length(all.clusters)
  
  if(K == 1)
  {
    edge.counts.matrix <- diag(1)
  }
  else
  {
    edge.counts.matrix <- Matrix(rep(0, K*K), c(K, K)) 
  }
  
  # loop through all clusters
  for(cluster.row in 1:K)
  {
    for(cluster.column in 1:K)
    {
      if(directed == FALSE)
      {
        edge.counts.matrix[cluster.row, cluster.column] <- 
          sum(adj[all.clusters[[cluster.row]], all.clusters[[cluster.column]]]) 
      }
      if(directed == TRUE)
      {
        edge.counts.matrix[cluster.row, cluster.column] <- 
          sum(adj[all.clusters[[cluster.row]], all.clusters[[cluster.column]]]) +
          sum(adj[all.clusters[[cluster.column]], all.clusters[[cluster.row]]]) 
      }
    }
  }
  return(edge.counts.matrix)
}
#CreateGlobalMatrixEdgeCountsBetweenClusters(all.clusters, directed)


#****************************************************
#'  Create global running total edge counts list
#'  
#' @param non.c.bar clusters that do not contain the anchors [list]
#' @param N number of particles [scalar]
#' @return edge counts [list of vectors]
CreateGlobalRunningTotalEdgeCountList <- function(non.c.bar, N)
{
  # IF no non.c.bar clusters
  if(length(non.c.bar) == 0)
  {
    edge.count.vector <- rep(0, 3)
  }
  
  # IF at least 1 non.c.bar cluster
  if(length(non.c.bar) > 0)
  {
    edge.count.vector <- rep(0, 3 + (2 * length(non.c.bar)))
  }
  
  previous.edges.list <- lapply(1:N, list)
  previous.edges.list <- lapply(previous.edges.list, function(x)
  {
    x <- edge.count.vector 
  })
  
  return(previous.edges.list)
}  
#CreateGlobalRunningTotalEdgeCountList(non.c.bar, N)

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
#'  Perform ancestor sampling
#' @param sigma Uniform permutation on closure of anchors (output from SamplePermutation) [vector]
#' @param s Anchor points [vector]
#' @param particle sequence of allocation decisions up to time t [vector]
#' @param log.previous.unnormalised.weight log of unnormalised weight at time t-1 [scalar]
#' @param log.gamma.hat.previous log intermediate target at time t-1 [scalar]
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
#' @param directed is network directed or not [boolean]
#' @param particle.index index of particle [scalar]
#' @param conditional.path the conditional path of allocation decisions (particle 1) [vector]
#' @return concatenated particle & ancestor sampling weight for the input particle [scalar]
AncestorSampling <- function(sigma, s, particle, log.previous.unnormalised.weight, 
                             log.gamma.hat.previous, all.clusters, non.c.bar, adj, 
                             tau1, tau2, t, n, alpha, beta1, beta2, directed, 
                             particle.index, conditional.path)
{
  # concatenate previous path of particle with future conditional path values
  concatenated.particle <- c(particle[1:(t-1)], conditional.path[t:n])
  
  # calculate log gamma hat for concatenated particle at final time t = n
  conc.log.gamma.hat <- LogImprovedIntermediateTarget(sigma, s, particle = concatenated.particle, 
                                                      all.clusters, non.c.bar, adj, tau1, tau2, 
                                                      t = n, n, alpha, beta1, beta2, directed, 
                                                      particle.index)
  
  # calculate ancestor sampling weights
  log.a.s.weight <- log.previous.unnormalised.weight + 
                    conc.log.gamma.hat$log.imp.int.target - 
                    log.gamma.hat.previous
  
  return(list("log.a.s.weight" = log.a.s.weight,
              "concatenated.particle" = concatenated.particle))
}


#****************************************************
#'  Particle Gibbs Split-Merge algorithm ("Algorithm 3")
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
#' @param directed is network directed or not? [boolean]
#' @param as.probability probability of ancester sampling step at each iteration [scalar]
#' @return Updated c.bar [list of vectors]
ParticleGibbsSplitMerge <- function(all.clusters, adj, s, s.bar, c.bar, non.c.bar, N, 
                                    resampling.threshold, alpha, beta1, beta2, directed,
                                    as.probability)
{
  # uniform permutation on elements of s.bar - with anchors 1st and 2nd
  sigma <- SamplePermutation(s, s.bar)
  n <- length(sigma)

  # define particle matrix & fix 1st particle of each generation to conditional path
  particles <- as.particles <- matrix(rep(0, N*n), c(N, n))
  particles[1,] <- conditional.path <- MapClustersToAllocations(sigma, c.bar) 
  particles[,1] <- rep(1, N) # column 1: 1st decision is (#1) initialise for all particles
  skip.particle.indices <- NULL # only need to calculate weights for single merge particle

  # define weights and previous log gamma hats at t=1
  log.un.weights <- rep(log(1), N)
  log.norm.weights <- rep(log(1/N), N) # 1st log norm weight is log(1/N) for all particles
  log.gamma.hat.previous <- rep(0, N)
  as.flag <- FALSE; 
  
  # set up global variables: for running totals of edge counts and log gamma at t=2
  global.running.total.edge.counts <<- CreateGlobalRunningTotalEdgeCountList(non.c.bar, N)
  global.running.total.max.counts <<- CreateGlobalRunningTotalEdgeCountList(non.c.bar, N)
  global.log.gamma_2 <<- rep(0, N)
  
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
    
    # skip weights calculation for any additional merge particles 
    particle.indices <- which(1:N %!in% skip.particle.indices == TRUE)
    
    # conditions required for ancestor sampling 
    if(t > 3 && runif(1) <= as.probability && conditional.path[2] != 2 && 
       num.split.particles > 1)
    {
      as.flag <- TRUE
      as.weights <- rep(0, N)
      as.particles <- array(rep(0, N * n), c(N, n))
    }
    
    for(p in particle.indices) # particle iterations
    {
      # ancestor sampling only if conditional path is split particle
      if(as.flag == TRUE && p %in% split.particle.indices)
      {
        as.output <- AncestorSampling(sigma, s, particle = particles[p, 1:t], 
                                      log.previous.unnormalised.weight = log.un.weights[p], 
                                      log.gamma.hat.previous = log.gamma.hat.previous[p], 
                                      all.clusters, non.c.bar, adj, tau1, tau2, t, n, alpha, 
                                      beta1, beta2, directed, particle.index = p, 
                                      conditional.path)
        as.weights[p] <- as.output$log.a.s.weight
        as.particles[p, ] <- as.output$concatenated.particle
      }
      
      # calculate proposal and weights
      weights.output <- LogUnnormalisedWeight(sigma, s, particle = particles[p, 1:t], 
                                              log.previous.unnormalised.weight = log.un.weights[p], 
                                              all.clusters, non.c.bar, adj, tau1, tau2, 
                                              t, n, alpha, beta1, beta2, directed,
                                              particle.index = p, log.gamma.hat.previous[p])
      
      # proposal allocations: don't change conditional path
      if(p >= 2)
      {
        particles[p, t] <- weights.output$proposal.allocation
      }
      
      # update weights and previous log gamma hat
      log.un.weights[p] <- weights.output$log.unnormalised.weight
      log.gamma.hat.previous[p] <- weights.output$log.gamma.hat 
    }
    
    # only calculate weights for 1 merge particle & check enough splits for A.S.
    if(t == 2)
    {
      # need at least 2 split particles for ancestor sampling
      split.particle.indices <- which(particles[,2] == 4)
      num.split.particles <- length(split.particle.indices)
      
      # only need to calculate weights for single merge particle
      merge.particle.indices <- which(particles[,2] == 2)
        
      if(length(merge.particle.indices) > 1)
      {
        first.merge.particle.index <- merge.particle.indices[1] 
        skip.particle.indices <- merge.particle.indices[-1]
      }
    }
    
    # select from ancestor sampling weights
    if(as.flag == TRUE)
    {
      # remove redundant zeros
      redundant.as.weight.indices <- which(as.weights == 0)
      non.redundant.as.weight.indices <- which(1:N %!in% redundant.as.weight.indices)
      as.weights <- as.weights[-redundant.as.weight.indices]
      
      # normalise weights and select an ancestral path
      log.norm.as.weights <- sapply(as.weights, function(x){x - logSumExp(as.weights)})
      select.path <- 
        which(as.double(rmultinom(n = 1, size = 1, prob = exp(log.norm.as.weights))) == 1)
      particles[1,] <- as.particles[non.redundant.as.weight.indices[select.path], ]
    }
    as.flag <- FALSE # reset as.flag at each iteration
  }
  
  # set weights for remaining merge particles  
  if(length(merge.particle.indices) > 1)
  {
    log.un.weights[skip.particle.indices] <- log.un.weights[first.merge.particle.index]
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
#' @param directed is network directed or not [boolean]
#' @param as.probability probability of ancester sampling step at each iteration [scalar]
#' @return Updated clustering
SplitMerge <- function(s, all.clusters, adj, N, resampling.threshold, alpha, 
                       beta1, beta2, directed, as.probability)
{
  #s <- SelectAnchors(all.clusters)  # select anchors uniformly

  # calculate c.bar and s.bar
  closure <- CalculateClosureOfAnchors(s, all.clusters)
  c.bar <- closure$c.bar  
  s.bar <- closure$s.bar
  
  # calculate non.c.bar
  non.c.bar.cluster.indicators <- lapply(lapply(all.clusters, 
                                                function(x){s %in% x}), 
                                         function(y){any(y == TRUE)})
  non.c.bar <- all.clusters[which(non.c.bar.cluster.indicators == FALSE)]
  
  # calculate global variables
  global.non.c.bar.indices <<- which(non.c.bar.cluster.indicators == FALSE)
  global.counts.between.nodes.clusters <<- 
    CreateGlobalMatrixEdgeCountsBetweenAllNodesAndClusters(all.clusters, num.nodes = dim(adj)[1], 
                                                           directed)
  #global.matrix.between.clusters.counts <<- 
  #  CreateGlobalMatrixEdgeCountsBetweenClusters(all.clusters, directed)
  
  # update clustering
  updated.c.bar <- ParticleGibbsSplitMerge(all.clusters, adj, s, s.bar, c.bar, non.c.bar, N, 
                                           resampling.threshold, alpha, beta1, beta2, directed,
                                           as.probability)
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

