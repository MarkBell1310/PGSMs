//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


//****************************************************
//' Update nxK (node:cluster) matrix for undirected network
//' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param cluster.to index of cluster that node is moved to [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated nxK matrix [matrix]
//[[Rcpp::export]]
arma::mat UpdateNxKMatrixUndirected(arma::mat prev_nxK_mat, int node_index, int cluster_from, 
                                    int cluster_to, int num_nodes, arma::mat adj_mat)
{
  // Undirected: we have only 1 nxK matrix.
  
  // Need to change indices in C++
  cluster_from = cluster_from - 1;
  cluster_to = cluster_to - 1;
  node_index = node_index - 1;
  
  // For each row, update the cluster_from and cluster_to columns
  for(int i = 0; i < num_nodes; i++)
  {               
    // (1) update "from" column - for cluster that node has left
    prev_nxK_mat(i, cluster_from) = prev_nxK_mat(i, cluster_from) - adj_mat(i, node_index);
    
    // (2) update "to" column - for cluster that node has joined
    prev_nxK_mat(i, cluster_to) = prev_nxK_mat(i, cluster_to) + adj_mat(i, node_index);
  }
  
  return prev_nxK_mat;
}


//****************************************************
//' Update nxK "TO" matrix: "from nodes TO CLUSTERS" for directed network
//' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
//' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param cluster.to index of cluster that node is moved to [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated nxK matrix [matrix]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::export]]
arma::mat UpdateNxKMatrixDirectedTo(arma::mat prev_nxK_mat_to, int node_index, int cluster_from, 
                                    int cluster_to, int num_nodes, arma::mat adj_mat)
{
  // Need to change indices in C++
  cluster_from = cluster_from - 1;
  cluster_to = cluster_to - 1;
  node_index = node_index - 1;
  
  // For each row, update the cluster_from and cluster_to columns
  for(int node = 0; node < num_nodes; node++)
  {
    // (1) update "from" column - for cluster that node has left
    prev_nxK_mat_to(node, cluster_from) = prev_nxK_mat_to(node, cluster_from) - adj_mat(node, node_index);
    
    // (2) update "to" column  - for cluster that node has joined
    prev_nxK_mat_to(node, cluster_to) = prev_nxK_mat_to(node, cluster_to) + adj_mat(node, node_index);
  }
  
  return prev_nxK_mat_to;
}


//****************************************************
//' Update nxK "FROM" matrix: "FROM CLUSTERS to nodes" for directed network
//' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
//' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param cluster.to index of cluster that node is moved to [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated nxK matrix [matrix]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::export]]
arma::mat UpdateNxKMatrixDirectedFrom(arma::mat prev_nxK_mat_from, int node_index, int cluster_from,
                                      int cluster_to, int num_nodes, arma::mat adj_mat)
{
  // Need to change indices in C++
  cluster_from = cluster_from - 1;
  cluster_to = cluster_to - 1;
  node_index = node_index - 1;
  
  // For each row, update the cluster_from and cluster_to columns
  for(int node = 0; node < num_nodes; node++)
  {
    // (1) update "from" column - for cluster that node has left
    prev_nxK_mat_from(node, cluster_from) = prev_nxK_mat_from(node, cluster_from) - 
      adj_mat(node_index, node);
    
    // (2) update "to" column - for cluster that node has joined
    prev_nxK_mat_from(node, cluster_to) = prev_nxK_mat_from(node, cluster_to) + adj_mat(node_index, node);
  }
  
  return prev_nxK_mat_from;
}


//****************************************************
//' Update extended nxK (node:cluster) matrix for undirected network
//' @param prev.nxK.mat nxK matrix from previous Gibbs iteration [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated extended nx(K+1) matrix [matrix]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::export]]
arma::mat UpdateExtendedNxKMatrixUndirected(arma::mat prev_nxK_mat, int node_index, int cluster_from, 
                                            int num_nodes, arma::mat adj_mat, int K)
{
  // Undirected: only 1 nxK matrix.
  
  // Create nx(K+1) matrix 
  arma::mat new_nxK_mat = arma::zeros(num_nodes, K+1);   
  
  // Need to change indices in C++
  int cluster_to = K + 1 - 1;         // node moves to K+1th cluster
  cluster_from = cluster_from - 1;
  node_index = node_index - 1;

  // For each row, update the cluster_from and cluster_to columns
  for(int node = 0; node < num_nodes; node++)
  {
    // (1) update "from" column - for cluster that node has left
    new_nxK_mat(node, cluster_from) = prev_nxK_mat(node, cluster_from) - adj_mat(node, node_index);

    // (2) update "to" column - for cluster that node has joined
    new_nxK_mat(node, cluster_to) = adj_mat(node, node_index);
  }
  
  // IF there are other columns to update in addition to cluster_from and cluster_to
  if(K + 1 > 2)
  {
    // (3) update any remaining columns using previous matrix
    for(int column = 0; column < K; column++)
    {
      if(column != cluster_from) // cluster_to is K+1th column so no need to include here
      {
        for(int node = 0; node < num_nodes; node++)
        {
          new_nxK_mat(node, column) = prev_nxK_mat(node, column);
        }
      }
    }
  }
    
  return new_nxK_mat;
}


//****************************************************
//' Update extended nxK "TO" matrix: "from nodes TO CLUSTERS" for directed network
//' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
//' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated extended nx(K+1) matrix [matrix]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::export]]
arma::mat UpdateExtendedNxKMatrixDirectedTo(arma::mat prev_nxK_mat_to, int node_index,
                                            int cluster_from, int num_nodes,
                                            arma::mat adj_mat, int K)
{
  // Directed: this function updates 1 of 2 nx(K+1) matrices

  // Create nx(K+1) matrix
  arma::mat new_nxK_mat_to = arma::zeros(num_nodes, K+1);

  // Need to change indices in C++
  int cluster_to = K + 1 - 1;         // node moves to K+1th cluster
  cluster_from = cluster_from - 1;
  node_index = node_index - 1;

  // For each row, update the cluster_from and cluster_to columns
  for(int node = 0; node < num_nodes; node++)
  {
    // (1) update "from" column - for cluster that node has left
    new_nxK_mat_to(node, cluster_from) = prev_nxK_mat_to(node, cluster_from) - adj_mat(node, node_index);
    
    // (2) update "to" column - for cluster that node has joined
    new_nxK_mat_to(node, cluster_to) = adj_mat(node, node_index);
  }

  // IF there are other columns to update in addition to cluster_from and cluster_to
  if(K + 1 > 2)
  {
    // (3) update any remaining columns using previous matrix
    for(int column = 0; column < K; column++)
    {
      if(column != cluster_from) // cluster_to is K+1th column so no need to include here
      {
        for(int node = 0; node < num_nodes; node++)
        {
          new_nxK_mat_to(node, column) = prev_nxK_mat_to(node, column);
        }
      }
    }
  }

  return new_nxK_mat_to;
}

//****************************************************
//' Update extended nxK "FROM" matrix: "FROM CLUSTERS to nodes" for directed network
//' @param prev.nxK.mat.to previous nxK "from nodes TO CLUSTERS" matrix [matrix]
//' @param prev.nxK.mat.from nxK previous nxK "FROM CLUSTERS to nodes" matrix [matrix]
//' @param node.index index of node being moved [scalar]
//' @param cluster.from index of cluster that node is moved from [scalar]
//' @param num.nodes number of nodes in network [scalar]
//' @param adj adjacency matrix [matrix]
//' @return Updated extended nx(K+1) matrix [matrix]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::export]]
arma::mat UpdateExtendedNxKMatrixDirectedFrom(arma::mat prev_nxK_mat_from, int node_index,
                                              int cluster_from, int num_nodes,
                                              arma::mat adj_mat, int K)
{
  // Directed: this function updates 2 of 2 nx(K+1) matrices

  // Create nx(K+1) matrix
  arma::mat new_nxK_mat_from = arma::zeros(num_nodes, K+1);

  // Need to change indices in C++
  int cluster_to = K + 1 - 1;         // node moves to K+1th cluster
  cluster_from = cluster_from - 1;
  node_index = node_index - 1;

  // For each row, update the cluster_from and cluster_to columns
  for(int node = 0; node < num_nodes; node++)
  {
    // (1) update "from" column - for cluster that node has left
    new_nxK_mat_from(node, cluster_from) = prev_nxK_mat_from(node, cluster_from) - adj_mat(node_index, node);
    
    // (2) update "to" column - for cluster that node has joined
    new_nxK_mat_from(node, cluster_to) = adj_mat(node_index, node);
  }

  // IF there are other columns to update in addition to cluster_from and cluster_to
  if(K + 1 > 2)
  {
    // (3) update any remaining columns using previous matrix
    for(int column = 0; column < K; column++)
    {
      if(column != cluster_from) // cluster_to is K+1th column so no need to include here
      {
        for(int node = 0; node < num_nodes; node++)
        {
          new_nxK_mat_from(node, column) = prev_nxK_mat_from(node, column);
        }
      }
    }
  }

  return new_nxK_mat_from;
}


