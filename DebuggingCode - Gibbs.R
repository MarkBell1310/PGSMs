

all.clusters = start.clusters

num.nodes.in.clusters <- previous.matrices$num.nodes.in.clusters
KxK.edge.counts <- previous.matrices$KxK.edge.counts
KxK.max.counts <- previous.matrices$KxK.max.counts

  
(M <- with( dat, 
            
undir <- sparseMatrix(i= 1:K, j= 1:K, x=rssi, dims = c(K, K))

i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))      


# if(k = 1st cluster, call the nxK and KxK edge count matrices, KxK max count matrix and vector
#    of num nodes in each cluster from  prev iter of either PGSMs/Gibbs)


InitialKxKEdgeCountMatrix <- function()