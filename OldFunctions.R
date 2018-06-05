
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
  # -the next m = length(c.bar[[2]]) elements are counts with c.bar cluster 2 (if it exists)

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
#PreComputeEdgeCountsBetweenCbarAndNonCbar(adj, non.c.bar, sigma, directed


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
