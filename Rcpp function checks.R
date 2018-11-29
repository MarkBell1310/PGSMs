
## Setup for checks
prev_nxK_mat = as.matrix(previous.matrices$nxK.edge.counts)
prev_nxK_mat_to = as.matrix(previous.matrices$nxK.edge.counts$edge.counts.to)
prev_nxK_mat_from = as.matrix(previous.matrices$nxK.edge.counts$edge.counts.from)
#prev_nxK_mat = matrix(rep(1:20, 3), c(20, 3))
K = 5
node_index = 9
cluster_from = 1
cluster_to = K+1
num_nodes = num.nodes
adj_mat = as.matrix(adj)
all.clusters


## nxK matrix undirected 
z <- UpdateNxKMatrixUndirectedR(prev.nxK.mat = previous.matrices$nxK.edge.counts, 
                                node.index = node_index, cluster.from = cluster_from, 
                                cluster.to = cluster_to, num.nodes = num_nodes, adj)

z - UpdateNxKMatrixUndirected(prev_nxK_mat, node_index, cluster_from, 
                              cluster_to, num_nodes, adj_mat)

## nxK matrices directed 
z <- UpdateNxKMatricesDirectedR(prev.nxK.mat.to = previous.matrices$nxK.edge.counts$edge.counts.to,
                          prev.nxK.mat.from = previous.matrices$nxK.edge.counts$edge.counts.from,
                          node.index = node_index, cluster.from = cluster_from, 
                          cluster.to = cluster_to, num.nodes, adj)

UpdateNxKMatrixDirectedTo(prev_nxK_mat_to, node_index, cluster_from, cluster_to, num_nodes, adj_mat) -
z$edge.counts.to 

UpdateNxKMatrixDirectedFrom(prev_nxK_mat_from, node_index, cluster_from, cluster_to, num_nodes, adj_mat) -
z$edge.counts.from

## Extended nxK matrix undirected
z <- UpdateExtendedNxKMatrixUndirected(prev_nxK_mat, node_index, cluster_from, 
                                       num_nodes, adj_mat, K)

z - UpdateExtendedNxKMatrixUndirectedR(prev.nxK.mat = prev_nxK_mat, node.index = node_index, 
                                   cluster.from = cluster_from, num.nodes = num_nodes, 
                                   adj, K)



## Extended nxK matrices directed
z <- UpdateExtendedNxKMatricesDirectedR(prev.nxK.mat.to = prev_nxK_mat_to, 
                                        prev.nxK.mat.from = prev_nxK_mat_from, 
                                        node.index = node_index, cluster.from = cluster_from, 
                                        num.nodes = num_nodes, adj, K)

z$edge.counts.to - UpdateExtendedNxKMatrixDirectedTo(prev_nxK_mat_to, node_index, cluster_from, 
                                                     num_nodes, adj_mat, K)
z$edge.counts.from - UpdateExtendedNxKMatrixDirectedFrom(prev_nxK_mat_from, node_index, cluster_from, 
                                                         num_nodes, adj_mat, K)


targ_nxK_mat = matrix(rep(0, 60), c(20, 3))
targ_nxK_mat[,1] = prev_nxK_mat[,1] - adj[,node_index]
targ_nxK_mat[,2] = prev_nxK_mat[,2] + adj[,node_index]
targ_nxK_mat[,3] = prev_nxK_mat[,3]

for(i in 1:5)
{
  print(prev_nxK_mat[,i] - z[,i])
}
