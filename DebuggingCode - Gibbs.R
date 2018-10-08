
#****************************************************
#************  DEBUG: Gibbs sampler *************
#****************************************************

## edge count checks
lapply(temp.storage.list, function(x){x$KxK.edge.counts})
lapply(temp.storage.list, function(x){x$KxK.max.counts})
lapply(temp.storage.list, function(x){x$num.nodes.in.clusters})
lapply(temp.storage.list, function(x){x$log.int.target})
       
sapply(temp.storage.list, function(x){any(x$KxK.edge.counts < 0)})
sapply(temp.storage.list, function(x){any(x$KxK.max.counts < 0)})
sapply(temp.storage.list, function(x){any(x$KxK.edge.counts > x$KxK.max.counts)})
sapply(temp.storage.list, function(x){any(sum(x$num.nodes.in.clusters) != num.nodes)})




temp.storage.list[[5]]$KxK.edge.counts
previous.matrices$KxK.edge.counts



all.clusters = start.clusters

num.nodes.in.clusters <- previous.matrices$num.nodes.in.clusters
KxK.edge.counts <- previous.matrices$KxK.edge.counts
KxK.max.counts <- previous.matrices$KxK.max.counts
