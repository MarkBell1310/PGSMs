
#****************************************************
#*********** Gibbs sampler - functions **************
#****************************************************

#****************************************************
#' Extend KxK matrix for K+1th cluster - add extra row and column
#' @param prev.mat KxK matrix from previous Gibbs iteration [matrix]
#' @param K number of clusters [scalar]
#' @return Extended (K+1)x(K+1) matrix [matrix]
ExtendKxKMatrix <- function(prev.matrix, K)
{
  cbind(rbind(prev.matrix, rep(0,K)), rep(0, K+1))
}

#****************************************************
#' Remove relevant rows/columns from edge count matrices when cluster.from is removed
#' @param previous.matrices edge count and max count matrices from previous iteration 
#' - (either from Gibbs or PGSMs)
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param num.nodes number of nodes in the network [scalar]
#' @param directed whether network is directed or not [boolean]
#' @return Updated matrices for loss of cluster.from [list of vectors and matrices]
RemoveClusterRowsColumns <- function(previous.matrices, cluster.from, num.nodes, directed)
{
  k <- cluster.from 
  
  # Update number of nodes in clusters
  previous.matrices$num.nodes.in.clusters <- previous.matrices$num.nodes.in.clusters[-cluster.from]
  
  # Remove kth row & column of edge count and max count matrices
  previous.matrices$KxK.edge.counts <- previous.matrices$KxK.edge.counts[-k, -k] 
  previous.matrices$KxK.max.counts <- previous.matrices$KxK.max.counts[-k, -k] 
  
  # Remove kth column of nxK matrix
  if(directed == FALSE)
  {
    previous.matrices$nxK.edge.counts <- previous.matrices$nxK.edge.counts[, -k]
  }
  if(directed == TRUE)
  {
    previous.matrices$nxK.edge.counts$edge.counts.from <- 
      previous.matrices$nxK.edge.counts$edge.counts.from[, -k]
    previous.matrices$nxK.edge.counts$edge.counts.to <- 
      previous.matrices$nxK.edge.counts$edge.counts.to[, -k]
  }
  
  # If now left with only 1 cluster, maintain matrix forms rather than scalars
  if(is.null(dim(previous.matrices$KxK.edge.counts)))
  {
    previous.matrices$KxK.edge.counts <- 
      sparseMatrix(1, 1, x = previous.matrices$KxK.edge.counts)
    previous.matrices$KxK.max.counts <- 
      sparseMatrix(1, 1, x = previous.matrices$KxK.max.counts)
    previous.matrices$nxK.edge.counts <- 
      sparseMatrix(1:num.nodes, rep(1, num.nodes), x = previous.matrices$nxK.edge.counts)
  }
  
  return(previous.matrices)
}

#****************************************************
#' Update KxK (cluster:cluster) edge counts matrix for undirected network
#' @param prev.KxK.mat KxK matrix from previous Gibbs iteration [matrix]
#' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param K number of clusters [scalar]
#' @return Updated KxK matrix [matrix]
UpdateKxKEdgeCountsMatrixUndirected <- function(prev.KxK.mat, prev.nxK.mat, node.index,
                                                cluster.from, cluster.to, K)
{
  ## Undirected: update relevant parts of (symmetric) upper triangular matrix only
  ## Note: We always have at least 2 clusters, since if K=1, then only need to use this function 
  ## for moving nodes to K+1th cluster - so no checks needed for types 1 and 2 counts

  new.KxK.mat <- prev.KxK.mat
  
  # "type 1" counts - within cluster (matrix diagonals)
  new.KxK.mat[cluster.from, cluster.from] <- 
    prev.KxK.mat[cluster.from, cluster.from] - prev.nxK.mat[node.index, cluster.from]
  
  new.KxK.mat[cluster.to, cluster.to] <- 
    prev.KxK.mat[cluster.to, cluster.to] + prev.nxK.mat[node.index, cluster.to]
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  # (the statements below ensure that the correct elements of upper triangular matrix
  #  are used and updated)
  if(cluster.from < cluster.to) 
  {
    new.KxK.mat[cluster.from, cluster.to] <- prev.KxK.mat[cluster.from, cluster.to] + 
      prev.nxK.mat[node.index, cluster.from] - prev.nxK.mat[node.index, cluster.to]
  }
  if(cluster.from > cluster.to)
  {
    new.KxK.mat[cluster.to, cluster.from] <- prev.KxK.mat[cluster.to, cluster.from] + 
      prev.nxK.mat[node.index, cluster.from] - prev.nxK.mat[node.index, cluster.to]
  }

  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 2)
  {
    for(k in (1:K)[-c(cluster.from, cluster.to)])
    {
      if(cluster.from < k) # as before
      {
        new.KxK.mat[cluster.from, k] <- prev.KxK.mat[cluster.from, k] - prev.nxK.mat[node.index, k]
      }
      if(cluster.from > k) # new
      {
        new.KxK.mat[k, cluster.from] <- prev.KxK.mat[k, cluster.from] - prev.nxK.mat[node.index, k]
      }
      
      if(cluster.to < k) # as before 
      {
        new.KxK.mat[cluster.to, k] <- prev.KxK.mat[cluster.to, k] + prev.nxK.mat[node.index, k]
      }
      if(cluster.to > k) # new
      {
        new.KxK.mat[k, cluster.to] <- prev.KxK.mat[k, cluster.to] + prev.nxK.mat[node.index, k]
      }
    }
  }
  
  # only need upper triangular part of matrix
  #new.KxK.mat[lower.tri(new.KxK.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.KxK.mat)
}

#****************************************************
#' Update KxK (cluster:cluster) edge counts matrix for directed network
#' @param prev.KxK.mat KxK matrix from previous Gibbs iteration [matrix]
#' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
#' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param K number of clusters [scalar]
#' @return Updated KxK matrix [matrix]
UpdateKxKEdgeCountsMatrixDirected <- function(prev.KxK.mat, prev.nxK.mat.to, prev.nxK.mat.from,
                                              node.index, cluster.from, cluster.to, K)
{
  ## Directed: need to update relevant parts of whole (asymmetric) matrix 
  ## Note: We always have at least 2 clusters, since if K=1, then only need to use this function 
  ## for moving nodes to K+1th cluster - so no checks needed for types 1 and 2 counts
  
  new.KxK.mat <- prev.KxK.mat
  
  # "type 1" counts - subtract both "to cluster" and "from cluster" edges
  new.KxK.mat[cluster.from, cluster.from] <- prev.KxK.mat[cluster.from, cluster.from] - 
    prev.nxK.mat.to[node.index, cluster.from] - prev.nxK.mat.from[node.index, cluster.from]
  
  new.KxK.mat[cluster.to, cluster.to] <- prev.KxK.mat[cluster.to, cluster.to] + 
    prev.nxK.mat.to[node.index, cluster.to] + prev.nxK.mat.from[node.index, cluster.to]
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  new.KxK.mat[cluster.from, cluster.to] <- prev.KxK.mat[cluster.from, cluster.to] -
    prev.nxK.mat.to[node.index, cluster.to] + prev.nxK.mat.from[node.index, cluster.from]
  
  new.KxK.mat[cluster.to, cluster.from] <- prev.KxK.mat[cluster.to, cluster.from] -
    prev.nxK.mat.from[node.index, cluster.to] + prev.nxK.mat.to[node.index, cluster.from]
  
  # "type 3" counts - between (Cf, others) and (Ct, others) - update both ways
  if(K > 2)
  {
    for(k in (1:K)[-c(cluster.from, cluster.to)])
    {
      new.KxK.mat[cluster.from, k] <- prev.KxK.mat[cluster.from, k] - prev.nxK.mat.to[node.index, k]
      new.KxK.mat[k, cluster.from] <- prev.KxK.mat[k, cluster.from] - prev.nxK.mat.from[node.index, k]
      
      new.KxK.mat[cluster.to, k] <- prev.KxK.mat[cluster.to, k] + prev.nxK.mat.to[node.index, k]
      new.KxK.mat[k, cluster.to] <- prev.KxK.mat[k, cluster.to] + prev.nxK.mat.from[node.index, k]
    }
  }
  return(new.KxK.mat)
}

#****************************************************
#' Calculate KxK matrix of maximum edge counts for undirected network 
#' @param prev.max.counts.mat KxK matrix of max counts from previous Gibbs iteration [matrix]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @return Updated KxK max counts matrix [matrix]
UpdateKxKMaxCountsMatrixUndirected <- function(prev.max.counts.mat, cluster.from, cluster.to, 
                                               K, all.clusters)
{
  ## Undirected: so update relevant parts of (symmetric) upper triangular matrix only
  #              -only rows/columns corresponding to cluster.from and cluster.to 
  new.max.counts.mat <- prev.max.counts.mat
  
  # Calculate num nodes in cluster.from, cluster.to WITHOUT the node that is moved
  # -note that all.clusters is the current clustering, so node has not moved yet 
  #   and so needs to be removed from cluster.from to calculate num nodes
  num.nodes.in.clusters <- sapply(all.clusters, length)
  num.nodes.cluster.from <- num.nodes.in.clusters[cluster.from] - 1
  num.nodes.cluster.to <- num.nodes.in.clusters[cluster.to] 
  
  # "type 1" counts - within cluster (matrix diagonals)
  new.max.counts.mat[cluster.from, cluster.from] <- 
    prev.max.counts.mat[cluster.from, cluster.from] -  num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.to] <- 
    prev.max.counts.mat[cluster.to, cluster.to] + num.nodes.cluster.to
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  # (the statements below ensure that the correct elements of upper triangular matrix
  #  are used and updated)
  if(cluster.from < cluster.to) 
  {
    new.max.counts.mat[cluster.from, cluster.to] <- 
      prev.max.counts.mat[cluster.from, cluster.to] - num.nodes.cluster.to + num.nodes.cluster.from
  }
  if(cluster.from > cluster.to)
  {
    new.max.counts.mat[cluster.to, cluster.from] <- 
      prev.max.counts.mat[cluster.to, cluster.from] - num.nodes.cluster.to + num.nodes.cluster.from
  }
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 2)
  {
    for(k in (1:K)[-c(cluster.from, cluster.to)])
    {
      if(cluster.from < k) 
      {
        new.max.counts.mat[cluster.from, k] <- prev.max.counts.mat[cluster.from, k] - 
          num.nodes.in.clusters[k]
      }
      if(cluster.from > k) 
      {
        new.max.counts.mat[k, cluster.from] <- prev.max.counts.mat[k, cluster.from] - 
          num.nodes.in.clusters[k]
      }

      if(cluster.to < k) 
      {
        new.max.counts.mat[cluster.to, k] <- prev.max.counts.mat[cluster.to, k] + 
          num.nodes.in.clusters[k]
      }
      if(cluster.to > k)
      {
        new.max.counts.mat[k, cluster.to] <- prev.max.counts.mat[k, cluster.to] + 
          num.nodes.in.clusters[k]
      }
    }
  }
  
  # only need upper triangular part of matrix
  #new.max.counts.mat[lower.tri(new.max.counts.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.max.counts.mat)
}

#****************************************************
#' Calculate KxK matrix of maximum edge counts for directed network 
#' @param prev.max.counts.mat KxK matrix of max counts from previous Gibbs iteration [matrix]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @return Updated KxK max counts matrix [matrix]
UpdateKxKMaxCountsMatrixDirected <- function(prev.max.counts.mat, cluster.from, cluster.to, 
                                               K, all.clusters)
{
  ## Directed: need to update relevant parts of whole (asymmetric) matrix 
  #            -only rows/columns corresponding to cluster.from and cluster.to 
  new.max.counts.mat <- prev.max.counts.mat
  
  # Calculate num nodes in cluster.from, cluster.to WITHOUT the node that is moved
  # -note that all.clusters is the current clustering, so node has not moved yet 
  #   and so needs to be removed from cluster.from to calculate num nodes
  num.nodes.in.clusters <- sapply(all.clusters, length)
  num.nodes.cluster.from <- num.nodes.in.clusters[cluster.from] - 1
  num.nodes.cluster.to <- num.nodes.in.clusters[cluster.to] 
  
  # "type 1" counts - within cluster (matrix diagonals)
  new.max.counts.mat[cluster.from, cluster.from] <- 
    prev.max.counts.mat[cluster.from, cluster.from] - 2 * num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.to] <- 
    prev.max.counts.mat[cluster.to, cluster.to] + 2 * num.nodes.cluster.to
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  new.max.counts.mat[cluster.from, cluster.to] <- 
    prev.max.counts.mat[cluster.from, cluster.to] - num.nodes.cluster.to + num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.from] <- new.max.counts.mat[cluster.from, cluster.to]
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 2)
  {
    for(k in (1:K)[-c(cluster.from, cluster.to)])
    {
      new.max.counts.mat[cluster.from, k] <- prev.max.counts.mat[cluster.from, k] - 
        num.nodes.in.clusters[k]
      new.max.counts.mat[k, cluster.from] <- new.max.counts.mat[cluster.from, k] 
      
      new.max.counts.mat[cluster.to, k] <- prev.max.counts.mat[cluster.to, k] + 
        num.nodes.in.clusters[k]
      new.max.counts.mat[k, cluster.to] <- new.max.counts.mat[cluster.to, k] 
    }
  }
  return(new.max.counts.mat)
}


#****************************************************
#' Update nxK (node:cluster) matrix for undirected network
#' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @return Updated nxK matrix [matrix]
UpdateNxKMatrixUndirected <- function(prev.nxK.mat, node.index, cluster.from, cluster.to, 
                                      num.nodes, adj)
{
  ## Undirected: only 1 nxK matrix
  new.nxK.mat <- prev.nxK.mat
  
  for(node in 1:num.nodes)
  {
    # (1) update "from" column - for cluster that node has left
    new.nxK.mat[node, cluster.from] <- prev.nxK.mat[node, cluster.from] - adj[node, node.index]
    
    # (2) update "to" column - for cluster that node has joined
    new.nxK.mat[node, cluster.to] <- prev.nxK.mat[node, cluster.to] + adj[node, node.index]
  }
  return(new.nxK.mat)
}

#****************************************************
#' Update nxK (node:cluster) matrices for directed network
#' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
#' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @return Updated nxK matrix [matrix]
UpdateNxKMatricesDirected <- function(prev.nxK.mat.to, prev.nxK.mat.from, node.index, 
                                      cluster.from, cluster.to, num.nodes, adj)
{
  ## Directed: 2 nxK matrices
  new.nxK.mat.to <- prev.nxK.mat.to
  new.nxK.mat.from <- prev.nxK.mat.from
  
  for(node in 1:num.nodes)
  {
    # (1) update "from" column of both matrices - for cluster that node has left
    new.nxK.mat.to[node, cluster.from] <- 
      prev.nxK.mat.to[node, cluster.from] - adj[node, node.index]
    new.nxK.mat.from[node, cluster.from] <- 
      prev.nxK.mat.from[node, cluster.from] - adj[node.index, node]
    
    # (2) update "to" column of both matrices - for cluster that node has joined
    new.nxK.mat.to[node, cluster.to] <- 
      prev.nxK.mat.to[node, cluster.to] + adj[node, node.index]
    new.nxK.mat.from[node, cluster.to] <- 
      prev.nxK.mat.from[node, cluster.to] + adj[node.index, node]
  }
  return(list("edge.counts.to" = new.nxK.mat.to,
              "edge.counts.from" = new.nxK.mat.from))
}


#****************************************************
#' Update matrices for clusters 1,...,K
#' @param previous.matrices edge count and max count matrices from previous iteration 
#' - (either from Gibbs or PGSMs)
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @param directed whether network is directed or not [boolean]
#' @return Updated matrices for kth cluster [list of matrices]
UpdateMatrices <- function(previous.matrices, node.index, cluster.from, cluster.to, K, 
                           all.clusters, num.nodes, adj, directed)
{
  # Undirected: update KxK edge counts, KxK max counts and single nxK matrix
  if(directed == FALSE)
  {
    return(list(
      "KxK.edge.counts" = 
        UpdateKxKEdgeCountsMatrixUndirected(prev.KxK.mat = previous.matrices$KxK.edge.counts, 
                                            prev.nxK.mat = previous.matrices$nxK.edge.counts, 
                                            node.index, cluster.from, cluster.to, K),
      "KxK.max.counts" = 
        UpdateKxKMaxCountsMatrixUndirected(prev.max.counts.mat = previous.matrices$KxK.max.counts, 
                                           cluster.from, cluster.to, K, all.clusters),
      "nxK.edge.counts" = 
        UpdateNxKMatrixUndirected(prev.nxK.mat = previous.matrices$nxK.edge.counts, 
                                  node.index, cluster.from, cluster.to, num.nodes, adj)))
  }
  
  # Directed: update KxK edge counts, KxK max counts and both nxK matrices
  if(directed == TRUE)
  {
    return(list(
      "KxK.edge.counts" = 
        UpdateKxKEdgeCountsMatrixDirected(prev.KxK.mat = previous.matrices$KxK.edge.counts, 
                                          prev.nxK.mat.to = previous.matrices$nxK.edge.counts$edge.counts.to, 
                                          prev.nxK.mat.from = previous.matrices$nxK.edge.counts$edge.counts.from,
                                          node.index, cluster.from, cluster.to, K),
      "KxK.max.counts" = 
        UpdateKxKMaxCountsMatrixDirected(prev.max.counts.mat = previous.matrices$KxK.max.counts, 
                                         cluster.from, cluster.to, K, all.clusters),
      "nxK.edge.counts" = 
      UpdateNxKMatricesDirected(prev.nxK.mat.to = previous.matrices$nxK.edge.counts$edge.counts.to,
                                prev.nxK.mat.from = previous.matrices$nxK.edge.counts$edge.counts.from,
                                node.index, cluster.from, cluster.to, num.nodes)))
  }
}

#****************************************************
#' Update K+1th row/column & attach to KxK edge counts matrix for undirected network
#' @param prev.KxK.mat KxK matrix from previous Gibbs iteration [matrix]
#' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param K number of clusters [scalar]
#' @return Updated extended (K+1)x(K+1) edge counts matrix [matrix]
UpdateExtendedKxKEdgeCountsUndirected <- function(prev.KxK.mat, prev.nxK.mat, node.index,
                                                  cluster.from, K)
{
  # Extend matrix by adding K+1th row and column
  new.KxK.mat <- ExtendKxKMatrix(prev.KxK.mat, K)
  cluster.to <- K+1
  
  ## Update "cluster from" (Cf) and "K+1th" (Ct) rows & columns ##
  
  # "type 1" counts - within clusters Cf And Ct (matrix diagonals)
  new.KxK.mat[cluster.from, cluster.from] <- 
    prev.KxK.mat[cluster.from, cluster.from] - prev.nxK.mat[node.index, cluster.from]
  new.KxK.mat[cluster.to, cluster.to] <- 0
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  if(cluster.from < cluster.to) 
  {
    new.KxK.mat[cluster.from, cluster.to] <- prev.nxK.mat[node.index, cluster.from]
  }
  if(cluster.from > cluster.to)
  {
    new.KxK.mat[cluster.to, cluster.from] <- prev.nxK.mat[node.index, cluster.from]
  }
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 1)
  {
    for(k in (1:K)[-cluster.from])
    {
      if(cluster.from < k) # as before
      {
        new.KxK.mat[cluster.from, k] <- prev.KxK.mat[cluster.from, k] - prev.nxK.mat[node.index, k]
      }
      if(cluster.from > k) # new
      {
        new.KxK.mat[k, cluster.from] <- prev.KxK.mat[k, cluster.from] - prev.nxK.mat[node.index, k]
      }

      if(cluster.to < k) # as before 
      {
        new.KxK.mat[cluster.to, k] <- prev.nxK.mat[node.index, k]
      }
      if(cluster.to > k) # new
      {
        new.KxK.mat[k, cluster.to] <- prev.nxK.mat[node.index, k]
      }
    }
  }
  
  # only need upper triangular part of matrix
  #new.KxK.mat[lower.tri(new.KxK.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.KxK.mat)
}

#****************************************************
#' Update K+1th row/column & attach to KxK edge counts matrix for directed network
#' @param prev.KxK.mat KxK matrix from previous Gibbs iteration [matrix]
#' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
#' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param K number of clusters [scalar]
#' @return Updated extended (K+1)x(K+1) edge counts matrix [matrix]
UpdateExtendedKxKEdgeCountsDirected <- function(prev.KxK.mat, prev.nxK.mat.to, prev.nxK.mat.from,
                                                node.index, cluster.from, K)
{
  # Extend matrix by adding K+1th row and column
  new.KxK.mat <- ExtendKxKMatrix(prev.KxK.mat, K)
  cluster.to <- K+1
  
  ## Update "cluster from" (Cf) and "K+1th" (Ct) rows & columns ##
  
  # "type 1" counts - within clusters Cf And Ct (matrix diagonals)
  new.KxK.mat[cluster.from, cluster.from] <- prev.KxK.mat[cluster.from, cluster.from] - 
    prev.nxK.mat.to[node.index, cluster.from] - prev.nxK.mat.from[node.index, cluster.from]
  new.KxK.mat[cluster.to, cluster.to] <- 0
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  new.KxK.mat[cluster.from, cluster.to] <- prev.nxK.mat.from[node.index, cluster.from]
  new.KxK.mat[cluster.to, cluster.from] <- prev.nxK.mat.to[node.index, cluster.from]
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 1)
  {
    for(k in (1:K)[-cluster.from])
    {
      new.KxK.mat[cluster.from, k] <- prev.KxK.mat[cluster.from, k] - prev.nxK.mat.to[node.index, k]
      new.KxK.mat[k, cluster.from] <- prev.KxK.mat[k, cluster.from] - prev.nxK.mat.from[node.index, k]
      
      new.KxK.mat[cluster.to, k] <- prev.nxK.mat.to[node.index, k]
      new.KxK.mat[k, cluster.to] <- prev.nxK.mat.from[node.index, k]
    }
  }
  
  # only need upper triangular part of matrix
  #new.KxK.mat[lower.tri(new.KxK.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.KxK.mat)
}

#****************************************************
#' Update K+1th row/column & attach to KxK max counts matrix for undirected network
#' @param prev.max.counts.mat KxK matrix of max counts from previous Gibbs iteration [matrix]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @return Updated extended (K+1)x(K+1) max counts matrix [matrix]
UpdateExtendedKxKMaxCountsUndirected <- function(prev.max.counts.mat, cluster.from, K,
                                                 all.clusters)
{
  # Extend matrix by adding K+1th row and column
  new.max.counts.mat <- ExtendKxKMatrix(prev.max.counts.mat, K)
  cluster.to <- K+1
  
  # Calculate num nodes in cluster.from, cluster.to WITHOUT the node that is moved
  # -note that all.clusters is the current clustering, so node has not moved yet 
  #   and so needs to be removed from cluster.from to calculate num nodes
  num.nodes.in.clusters <- sapply(all.clusters, length)
  num.nodes.cluster.from <- num.nodes.in.clusters[cluster.from] - 1
  num.nodes.cluster.to <- 0
  
  # "type 1" counts - within cluster (matrix diagonals)
  new.max.counts.mat[cluster.from, cluster.from] <- 
    prev.max.counts.mat[cluster.from, cluster.from] - num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.to] <- 0
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  if(cluster.from < cluster.to) 
  {
    new.max.counts.mat[cluster.from, cluster.to] <- num.nodes.cluster.from
  }
  if(cluster.from > cluster.to)
  {
    new.max.counts.mat[cluster.to, cluster.from] <- num.nodes.cluster.from
  }
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 1)
  {
    for(k in (1:K)[-cluster.from])
    {
      if(cluster.from < k) 
      {
        new.max.counts.mat[cluster.from, k] <- prev.max.counts.mat[cluster.from, k] - 
          num.nodes.in.clusters[k]
      }
      if(cluster.from > k) 
      {
        new.max.counts.mat[k, cluster.from] <- prev.max.counts.mat[k, cluster.from] - 
          num.nodes.in.clusters[k]
      }
      
      if(cluster.to < k) 
      {
        new.max.counts.mat[cluster.to, k] <- num.nodes.in.clusters[k]
      }
      if(cluster.to > k) 
      {
        new.max.counts.mat[k, cluster.to] <- num.nodes.in.clusters[k]
      }
    }
  }

  # only need upper triangular part of matrix
  #new.max.counts.mat[lower.tri(new.max.counts.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.max.counts.mat)
}

#****************************************************
#' Update K+1th row/column & attach to KxK max counts matrix for directed network
#' @param prev.max.counts.mat KxK matrix of max counts from previous Gibbs iteration [matrix]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @return Updated extended (K+1)x(K+1) max counts matrix [matrix]
UpdateExtendedKxKMaxCountsDirected <- function(prev.max.counts.mat, cluster.from, K,
                                               all.clusters)
{
  # Extend matrix by adding K+1th row and column
  new.max.counts.mat <- ExtendKxKMatrix(prev.max.counts.mat, K)
  cluster.to <- K+1
  
  # Calculate num nodes in cluster.from, cluster.to WITHOUT the node that is moved
  # -note that all.clusters is the current clustering, so node has not moved yet 
  #   and so needs to be removed from cluster.from to calculate num nodes
  num.nodes.in.clusters <- sapply(all.clusters, length)
  num.nodes.cluster.from <- num.nodes.in.clusters[cluster.from] - 1
  num.nodes.cluster.to <- 0
  
  # "type 1" counts - within cluster (matrix diagonals)
  new.max.counts.mat[cluster.from, cluster.from] <- 
    prev.max.counts.mat[cluster.from, cluster.from] - 2 * num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.to] <- 0
  
  # "type 2" counts - between Cf and Ct (cluster.from and cluster.to)
  new.max.counts.mat[cluster.from, cluster.to] <- num.nodes.cluster.from
  new.max.counts.mat[cluster.to, cluster.from] <- num.nodes.cluster.from
  
  # "type 3" counts - between (Cf, others) and (Ct, others)
  if(K > 1)
  {
    for(k in (1:K)[-cluster.from])
    {
      new.max.counts.mat[cluster.from, k] <- prev.max.counts.mat[cluster.from, k] - 
        num.nodes.in.clusters[k]
      new.max.counts.mat[k, cluster.from] <- prev.max.counts.mat[k, cluster.from] - 
        num.nodes.in.clusters[k]
      
      new.max.counts.mat[cluster.to, k] <- num.nodes.in.clusters[k]
      new.max.counts.mat[k, cluster.to] <- num.nodes.in.clusters[k]
    }
  }
  
  # only need upper triangular part of matrix
  #new.max.counts.mat[lower.tri(new.max.counts.mat)] <- rep(0, K*(K+1)/2)
  
  return(new.max.counts.mat)
}

#****************************************************
#' Update extended nxK (node:cluster) matrix for undirected network
#' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @return Updated extended nx(K+1) matrix [matrix]
UpdateExtendedNxKMatrixUndirected <- function(prev.nxK.mat, node.index, cluster.from, 
                                              num.nodes, adj, K)
{
  ## Undirected: only 1 nxK matrix - Extend matrix by adding K+1th column
  new.nxK.mat <- cbind(prev.nxK.mat, rep(0, num.nodes))
  cluster.to <- K+1
  
  for(node in 1:num.nodes)
  {
    # (1) update "from" column - for cluster that node has left
    new.nxK.mat[node, cluster.from] <- prev.nxK.mat[node, cluster.from] - adj[node, node.index]
    
    # (2) update "to" column - for cluster that node has joined
    new.nxK.mat[node, cluster.to] <- adj[node, node.index]
  }
  return(new.nxK.mat)
}

#****************************************************
#' Update extended nxK (node:cluster) matrices for directed network
#' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
#' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @return Updated extended nx(K+1) matrix [matrix]
UpdateExtendedNxKMatricesDirected <- function(prev.nxK.mat.to, prev.nxK.mat.from, node.index, 
                                            cluster.from, num.nodes, adj, K)
{
  ## Directed: need 2 nxK matrices - Extend matrices by adding K+1th columns
  new.nxK.mat.to <- cbind(prev.nxK.mat.to, rep(0, num.nodes))
  new.nxK.mat.from <- cbind(prev.nxK.mat.from, rep(0, num.nodes))
  cluster.to <- K+1
  
  for(node in 1:num.nodes)
  {
    # (1) update "from" column of both matrices - for cluster that node has left
    new.nxK.mat.to[node, cluster.from] <- 
      prev.nxK.mat.to[node, cluster.from] - adj[node, node.index]
    new.nxK.mat.from[node, cluster.from] <- 
      prev.nxK.mat.from[node, cluster.from] - adj[node.index, node]
    
    # (2) update "to" column of both matrices - for cluster that node has joined
    new.nxK.mat.to[node, cluster.to] <- adj[node, node.index]
    new.nxK.mat.from[node, cluster.to] <- adj[node.index, node]
  }
  return(list("edge.counts.to" = new.nxK.mat.to,
              "edge.counts.from" = new.nxK.mat.from))
}


#****************************************************
#' Update extended matrices for K+1th cluster
#' @param previous.matrices edge count and max count matrices from previous iteration 
#' - (either from Gibbs or PGSMs)
#' @param node.index index of node being moved [scalar]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param K number of clusters [scalar]
#' @param all.clusters current clustering [list of vectors]
#' @param num.nodes number of nodes in network [scalar]
#' @param adj adjacency matrix [matrix]
#' @param directed whether network is directed or not [boolean]
#' @return Updated extended matrices for K+1th cluster [list of matrices]
UpdateExtendedMatrices <- function(previous.matrices, node.index, cluster.from, 
                                   K, all.clusters, num.nodes, adj, directed)
{
  # Undirected: update KxK edge counts, KxK max counts and single nxK matrix
  if(directed == FALSE)
  {
    return(list(
      "KxK.edge.counts" = 
        UpdateExtendedKxKEdgeCountsUndirected(prev.KxK.mat = previous.matrices$KxK.edge.counts, 
                                              prev.nxK.mat = previous.matrices$nxK.edge.counts,
                                              node.index, cluster.from, K),
      "KxK.max.counts" = 
        UpdateExtendedKxKMaxCountsUndirected(prev.max.counts.mat = previous.matrices$KxK.max.counts, 
                                             cluster.from, K, all.clusters),
      "nxK.edge.counts" = 
        UpdateExtendedNxKMatrixUndirected(prev.nxK.mat = previous.matrices$nxK.edge.counts,
                                          node.index, cluster.from, num.nodes, adj, K)))
  }
  
  # Directed: update KxK edge counts, KxK max counts and both nxK matrices
  if(directed == TRUE)
  {
    return(list(
      "KxK.edge.counts" = 
        UpdateExtendedKxKEdgeCountsDirected(prev.KxK.mat = previous.matrices$KxK.edge.counts, 
                                            prev.nxK.mat.to = previous.matrices$nxK.edge.counts$edge.counts.to, 
                                            prev.nxK.mat.from = previous.matrices$nxK.edge.counts$edge.counts.from,
                                            node.index, cluster.from, K),
      "KxK.max.counts" = 
        UpdateExtendedKxKMaxCountsDirected(prev.max.counts.mat = previous.matrices$KxK.max.counts, 
                                           cluster.from, K, all.clusters),
      "nxK.edge.counts" = 
        UpdateExtendedNxKMatricesDirected(prev.nxK.mat.to = previous.matrices$nxK.edge.counts$edge.counts.to,
                                          prev.nxK.mat.from = previous.matrices$nxK.edge.counts$edge.counts.from, 
                                          node.index, cluster.from, num.nodes, adj, K)))
  }
}

#****************************************************
#' Update number of nodes in each cluster - when node moves to EXISTING cluster
#' @param prev.num.nodes.in.clusters previous num nodes in each cluster [vector]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @return Updated number of nodes in each cluster [vector]
UpdateNumNodesInClusters <- function(prev.num.nodes.in.clusters, cluster.from, cluster.to)
{
  # adjust counts for node that has moved 
  new.num.nodes.in.clusters <- prev.num.nodes.in.clusters
  new.num.nodes.in.clusters[cluster.from] <- prev.num.nodes.in.clusters[cluster.from] - 1
  new.num.nodes.in.clusters[cluster.to] <- prev.num.nodes.in.clusters[cluster.to] + 1
  return(new.num.nodes.in.clusters)
}

#****************************************************
#' Update number of nodes in each cluster - when node moves to NEW cluster
#' @param prev.num.nodes.in.clusters previous num nodes in each cluster [vector]
#' @param cluster.from index of cluster that node is moved from [scalar]
#' @param cluster.to index of cluster that node is moved to [scalar]
#' @return Updated number of nodes in each cluster [vector]
UpdateNumNodesInClustersExtended <- function(prev.num.nodes.in.clusters, cluster.from, cluster.to)
{
  # adjust counts for node that has moved 
  new.num.nodes.in.clusters <- prev.num.nodes.in.clusters
  new.num.nodes.in.clusters[cluster.from] <- prev.num.nodes.in.clusters[cluster.from] - 1
  new.num.nodes.in.clusters <- c(new.num.nodes.in.clusters, 1)
  return(new.num.nodes.in.clusters)
}

#****************************************************
#' Log intermediate target for Gibbs sampler
#' @param KxK.edge.counts KxK edge counts matrix [matrix]
#' @param KxK.max.counts KxK maximum counts matrix [matrix]
#' @param num.nodes.in.clusters number of nodes in each cluster [vector]
#' @param alpha parameter for Dirichlet [scalar]
#' @param beta1 likelihood tuning parameter [scalar]
#' @param beta2 likelihood tuning parameter [scalar]
#' @param K number of clusters [scalar]
#' @param num.nodes number of nodes [scalar]
#' @param directed whether network is directed or not [boolean]
#' @return log intermediate target [scalar]
LogIntermediateTargetGibbs <- function(KxK.edge.counts, KxK.max.counts, num.nodes.in.clusters, 
                                       alpha, beta1, beta2, K, num.nodes, directed)
{
  # sum over K clusters of log(gamma(num.nodes.in.cluster + alpha))
  sum.over.clusters <- sum(sapply(num.nodes.in.clusters, function(x)
  {
    lgamma(x + alpha)
  }))
  
  # define edge and max counts
  if(directed == TRUE)
  {
    ## Directed: both upper and lower triangular
    edge.counts <- as.vector(KxK.edge.counts)
    max.counts <- as.vector(KxK.max.counts)
  }
  else
  {
    ## Undirected: upper triangular only
    edge.counts <- KxK.edge.counts[upper.tri(KxK.edge.counts, diag = TRUE)]
    max.counts <- KxK.max.counts[upper.tri(KxK.max.counts, diag = TRUE)]
  }
  
  # sum over blocks of log-likelihoods: log(f(x_(kl)|z))
  sum.log.likelihoods <- sum(sapply(Map(function(x, y)
  {
    lbeta(beta1 + x, y - x + beta2) - lbeta(beta1, beta2)
  }, 
  edge.counts, max.counts), function(x){x}))

  # log intermediate target (Uses the prior in (6) of Mcdaid et al. (2013))
   - lfactorial(K) + lgamma(alpha * K) + sum.over.clusters - K * lgamma(alpha) - 
    lgamma(num.nodes + (alpha * K)) + sum.log.likelihoods
  #- log(factorial(K)) + log(gamma(alpha * K)) + sum.over.clusters - K * log(gamma(alpha)) - 
  #  log(gamma(num.nodes + (alpha * K))) + sum.log.likelihoods
}

#****************************************************
#' Full Gibbs sweep
#' @param all.clusters Current clustering [list of vectors]
#' @param alpha parameter for Dirichlet [scalar]
#' @param beta1 likelihood tuning parameter [scalar]
#' @param beta2 likelihood tuning parameter [scalar]
#' @param num.nodes number of nodes [scalar]
#' @param previous.matrices edge count and max count matrices from previous iteration 
#' - (either from Gibbs or PGSMs)
#' @param iter iteration number [scalar]
#' @param directed whether network is directed or not [boolean]
#'# @param global.counts.between.nodes.clusters Edge counts between all nodes & clusters [matrix]
#' @return Updated clustering and previous matrices [list]
GibbsSweep <- function(all.clusters, alpha, beta1, beta2, num.nodes, previous.matrices,
                       iter, directed)
{
  # Perform deterministic scan but with a random order
  random.nodes.order <- sample(1:num.nodes)
  
  ## (Do a Gibbs iteration and try node in each of the K+1 clusters) ##
  for(node in random.nodes.order)
  {
    #node <- random.nodes.order[3]
    K <- length(all.clusters)
    cluster.from <- which(lapply(all.clusters, function(x){which(node %in% x)}) == 1)

    # if cluster.from only contains 1 node, only need to investigate K clusters
    if(length(all.clusters[[cluster.from]]) == 1)
    {
      temp.storage.list <- lapply(1:K, list)
    
    } else 
    {
      temp.storage.list <- lapply(1:(K+1), list)
    }
    
    ## (Try node in the kth cluster) ##
    for(k in 1:(K+1))      
    {
      cluster.to <- k
      
      # if node already in proposed cluster, use previous matrices and log gamma
      if(cluster.to == cluster.from)
      {
        temp.storage.list[[k]] <- list("num.nodes.in.clusters" = previous.matrices$num.nodes.in.clusters,
                                       "KxK.max.counts" = previous.matrices$KxK.max.counts,
                                       "KxK.edge.counts" = previous.matrices$KxK.edge.counts,
                                       "nxK.edge.counts" = previous.matrices$nxK.edge.counts,
                                       "log.int.target" = previous.matrices$log.int.target,
                                       "cluster.to" = cluster.from)
        next
      }
      
      ##*****************************************************
      ## Update matrices and log intermediate target - then store results temporarily
      
      # if node moves to existing cluster
      if(k < K+1)
      {
        updated.matrices <- UpdateMatrices(previous.matrices, node.index = node, cluster.from, 
                                           cluster.to, K, all.clusters, num.nodes, adj, directed)
        num.nodes.in.clusters <- 
          UpdateNumNodesInClusters(prev.num.nodes.in.clusters = previous.matrices$num.nodes.in.clusters, 
                                   cluster.from, cluster.to)
        
        # check if K decreases after node moves - need appropriate K for log.int.target
        if(num.nodes.in.clusters[cluster.from] == 0)
        {
          temp.num.clusters <- K-1
        
        } else 
        {
          temp.num.clusters <- K
        }
        
        log.int.target <-
          LogIntermediateTargetGibbs(KxK.edge.counts = updated.matrices$KxK.edge.counts,
                                     KxK.max.counts = updated.matrices$KxK.max.counts,
                                     num.nodes.in.clusters, alpha, beta1, beta2,
                                     K = temp.num.clusters, num.nodes, directed)
        
        # Store updated nxK, KxK matrices and number of nodes in each cluster temporarily
        temp.storage.list[[k]] <- list("num.nodes.in.clusters" = num.nodes.in.clusters,
                                       "KxK.max.counts" = updated.matrices$KxK.max.counts,
                                       "KxK.edge.counts" = updated.matrices$KxK.edge.counts,
                                       "nxK.edge.counts" = updated.matrices$nxK.edge.counts,
                                       "log.int.target" = log.int.target,
                                       "cluster.to" = cluster.to)
      }
      
      # if node moves to new cluster - dimension of matrices increases
      if(k == K+1 && length(all.clusters[[cluster.from]]) > 1)
      {
        updated.matrices <- UpdateExtendedMatrices(previous.matrices, node.index = node, 
                                                   cluster.from, K, all.clusters, 
                                                   num.nodes, adj, directed)
        num.nodes.in.clusters <- 
          UpdateNumNodesInClustersExtended(prev.num.nodes.in.clusters = previous.matrices$num.nodes.in.clusters, 
                                           cluster.from, cluster.to)
        log.int.target <- 
          LogIntermediateTargetGibbs(KxK.edge.counts = updated.matrices$KxK.edge.counts, 
                                     KxK.max.counts = updated.matrices$KxK.max.counts, 
                                     num.nodes.in.clusters, alpha, beta1, beta2, 
                                     K = K+1, num.nodes, directed)
        
        # Store updated nxK, KxK matrices and number of nodes in each cluster temporarily
        temp.storage.list[[k]] <- list("num.nodes.in.clusters" = num.nodes.in.clusters,
                                       "KxK.max.counts" = updated.matrices$KxK.max.counts,
                                       "KxK.edge.counts" = updated.matrices$KxK.edge.counts,
                                       "nxK.edge.counts" = updated.matrices$nxK.edge.counts,
                                       "log.int.target" = log.int.target,
                                       "cluster.to" = cluster.to)
      }
    }
    
    ##*****************************************************
    ## Determine node location & store relevant matrices for next Gibbs iteration 
    
    # Normalise the log intermediate targets ("log gammas") 
    log.gammas <- sapply(temp.storage.list, function(x){x$log.int.target})
    log.norm.gammas <- sapply(log.gammas, function(x){x - logSumExp(log.gammas)})
    
    # Sample from these to determine which cluster the node should be moved to
    multinomial.sample <- as.double(rmultinom(n = 1, size = 1, prob = exp(log.norm.gammas)))
    
    # Store the edge count matrices that represent where the node was moved to 
    # - these represent the updated edge counts for the iteration involving the next node
    previous.matrices <- temp.storage.list[[which(multinomial.sample == 1)]]
    cluster.move.to <- previous.matrices$cluster.to
    
    ##*****************************************************
    ## Update clustering
    
    if(previous.matrices$cluster.to != cluster.from)
    {
      # Remove moving node from old cluster 
      all.clusters[[cluster.from]] <- 
        all.clusters[[cluster.from]][-which(all.clusters[[cluster.from]] == node)]
      
      # Move node to new cluster 
      if(cluster.move.to < K+1)
      {
        all.clusters[[cluster.move.to]] <- c(all.clusters[[cluster.move.to]], node)
      }
      if(cluster.move.to == K+1)
      {
        all.clusters[[cluster.move.to]] <- node
      }
      
      # Remove cluster.from if contains no nodes & remove relevant rows/cols of previous.matrices
      if(length(all.clusters[[cluster.from]]) == 0)
      {
        all.clusters[[cluster.from]] <- NULL
        previous.matrices <- RemoveClusterRowsColumns(previous.matrices, cluster.from, 
                                                      num.nodes, directed)
      }
    }
    #all.clusters
  }
  
  return(list("previous.matrices" = previous.matrices,
              "all.clusters" = all.clusters))
}


