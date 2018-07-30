
#****************************************************
#*********** Initialisation functions ***************
#****************************************************
#
# Functions in this script set up and calculate the inital vectors and matrices 
# that will be required (and updated) by both the PGSMs and the Gibbs sampler
#
# Note: All are used only once at start of algorithm and so are not optimised for speed.
#
#**************************************************** 

#****************************************************
#' Initial nxK (node:cluster) matrices for undirected AND directed network
#' @param all.clusters all current clusters [list of vectors]
#' @param num.nodes number of nodes [scalar]
#' @param adj adjacency matrix [matrix]
#' @param directed whether netowrk is directed or not [boolean]
#' @return edge counts [matrix: dim = num.nodes x K]
InitialNxKMatrices <- function(all.clusters, num.nodes, adj, directed)
{
  K <- length(all.clusters)
  
  # Undirected: only 1 nxK matrix required
  if(directed == FALSE)
  {
    edge.counts.matrix <- Matrix(rep(0, num.nodes*K), c(num.nodes, K), sparse = TRUE)
    
    # loop through all nodes and all clusters
    for(node in 1:num.nodes)
    {
      for(cluster in 1:K)
      {
        edge.counts.matrix[node, cluster] <- sum(adj[node, all.clusters[[cluster]]])
      }
    }
    return(edge.counts.matrix)
  }
  
  # Directed: 2 nxK matrices required: "to" and from" 
  # "to" = "from nodes TO CLUSTERS" 
  # "from" = "FROM CLUSTERS to nodes" 
  if(directed == TRUE)
  {
    edge.counts.to <- edge.counts.from <- Matrix(rep(0, num.nodes*K), c(num.nodes, K), 
                                                 sparse = TRUE)
    # loop through all nodes and all clusters
    for(node in 1:num.nodes)
    {
      for(cluster in 1:K)
      {
        edge.counts.to[node, cluster] <- sum(adj[node, all.clusters[[cluster]]])
        edge.counts.from[node, cluster] <- sum(adj[all.clusters[[cluster]]], node)
      }
    }
    return(list("edge.counts.to" = edge.counts.to,
                "edge.counts.from" = edge.counts.from))
  }
}
#InitialNxKMatrices(all.clusters, num.nodes, adj, directed)


#****************************************************
#' Initial KxK (cluster:cluster) edge count matrix for undirected AND directed network
#' @param all.clusters all current clusters [list of vectors]
#' @param adj adjacency matrix [matrix]
#' @param directed whether netowrk is directed or not [boolean]
#' @return edge counts [matrix: dim = K x K]
InitialKxKEdgeCountMatrix <- function(all.clusters, adj, directed)
{
  K <- length(all.clusters)
  
  # define sparse matrix
  if(K == 1)
  {
    edge.counts.matrix <- Matrix(diag(K), sparse = TRUE)
  }
  else
  {
    edge.counts.matrix <- Matrix(rep(0, K*K), c(K, K), sparse = TRUE)
  }
  
  # calculate edge counts by looping through all clusters
  for(cluster.row in 1:K)
  {
    for(cluster.column in 1:K)
    {
      edge.counts.matrix[cluster.row, cluster.column] <- 
        sum(adj[all.clusters[[cluster.row]], all.clusters[[cluster.column]]]) 
    }
  }
  
  if(directed == FALSE)
  {
    # only need upper triangular part of matrix
    edge.counts.matrix[lower.tri(edge.counts.matrix)] <- rep(0, K*(K+1)/2)
  }
  
  return(edge.counts.matrix)
}
#InitialKxKEdgeCountMatrix(all.clusters, adj)

#****************************************************
#' Initial number of nodes in each cluster
#' @param all.clusters all current clusters [list of vectors]
#' @return Initial number of nodes in each cluster [vector]
InitialNumNodesInClusters <- function(all.clusters)
{
  sapply(all.clusters, length)
}
#InitialNumNodesInClusters(all.clusters)

#****************************************************
#' Initial KxK (cluster:cluster) maximum count matrix for undirected AND directed network
#' @param all.clusters all current clusters [list of vectors]
#' @param num.nodes.in.clusters number of nodes in each cluster [vector]
#' @param directed whether netowrk is directed or not [boolean]
#' @return maximum counts [matrix: dim = K x K]
InitialKxKMaxCountMatrix <- function(all.clusters, num.nodes.in.clusters, directed)
{
  K <- length(all.clusters)
  
  # define sparse matrix
  if(K == 1)
  {
    max.counts.matrix <- Matrix(diag(K), sparse = TRUE)
  }
  else
  {
    max.counts.matrix <- Matrix(rep(0, K*K), c(K, K), sparse = TRUE)
  }
  
  # loop through all clusters
  for(cluster.row in 1:K)
  {
    for(cluster.column in 1:K)
    {
      # off-diagonal counts are same for undirected/directed network
      if(cluster.row != cluster.column)
      {
        max.counts.matrix[cluster.row, cluster.column] <- 
          num.nodes.in.clusters[cluster.row] * num.nodes.in.clusters[cluster.column]
      }
      
      # diagonal counts are different for undirected/directed network
      else
      {
        if(directed == FALSE)
        {
          max.counts.matrix[cluster.row, cluster.column] <- 
            0.5 * num.nodes.in.clusters[cluster.row] * (num.nodes.in.clusters[cluster.row] - 1)
        }
        if(directed == TRUE)
        {
          max.counts.matrix[cluster.row, cluster.column] <- 
            num.nodes.in.clusters[cluster.row] * (num.nodes.in.clusters[cluster.row] - 1)
        }
      }
    }
  }
  
  if(directed == FALSE)
  {
    # only need upper triangular part of matrix
    max.counts.matrix[lower.tri(max.counts.matrix)] <- rep(0, K*(K+1)/2)
  }
  
  return(max.counts.matrix)
}
#InitialKxKMaxCountMatrix(all.clusters, num.nodes.in.clusters, directed)

#****************************************************
#' Initial KxK (cluster:cluster) maximum count matrix for undirected AND directed network
#' @param all.clusters all current clusters [list of vectors]
#' @param num.nodes number of nodes [scalar]
#' @param adj adjacency matrix [matrix]
#' @param directed whether network is directed or not [boolean]
#' @return maximum counts [matrix: dim = K x K]
InitialSetupList <- function(all.clusters, num.nodes, adj, directed)
{
  num.nodes.in.clusters <- InitialNumNodesInClusters(all.clusters)
  KxK.max.counts <- InitialKxKMaxCountMatrix(all.clusters, num.nodes.in.clusters, directed)
  KxK.edge.counts <- InitialKxKEdgeCountMatrix(all.clusters, adj, directed)
  nxK.edge.counts <- InitialNxKMatrices(all.clusters, num.nodes, adj, directed)
  
  return(list("num.nodes.in.clusters" = num.nodes.in.clusters,
              "KxK.max.counts" = KxK.max.counts,
              "KxK.edge.counts" = KxK.edge.counts,
              "nxK.edge.counts" = nxK.edge.counts))
}
#InitialSetupList(all.clusters, num.nodes, adj, directed)
