

previous.matrices$KxK.edge.counts
previous.matrices$KxK.max.counts
previous.matrices$nxK.edge.counts$edge.counts.to
previous.matrices$nxK.edge.counts$edge.counts.from

## KxK edge count checks
diags <- rep(0, 5)
for(i in 1:5)
{
  diags[i] <- sum(adj[all.clusters[[i]], all.clusters[[i]]])
}
sum(adj[all.clusters[[1]], all.clusters[[2]]])
sum(adj[all.clusters[[2]], all.clusters[[6]]])

sum(adj[all.clusters[[1]], all.clusters[[6]]])
sum(adj[all.clusters[[6]], all.clusters[[1]]])

sum(adj[all.clusters[[3]], all.clusters[[6]]])
sum(adj[all.clusters[[4]], all.clusters[[6]]])

sum(adj[all.clusters[[5]], all.clusters[[2]]])
sum(adj[all.clusters[[3]], all.clusters[[2]]])

## KxK max count checks
all.clusters
sapply(all.clusters, function(x) {MaxNodesWithinCluster(length(x), directed = TRUE)})
length(all.clusters[[x]]) * length(all.clusters[[y]])

## n x K matrices checks
vec1 <- vec2 <- vec3 <- vec4 <- vec5 <- vec6 <- rep(0, 20);  
for(i in 1:20)
{
  vec1[i] <- sum(adj[i, all.clusters[[1]]])
  vec2[i] <- sum(adj[i, all.clusters[[2]]])
  vec3[i] <- sum(adj[i, all.clusters[[3]]])
  vec4[i] <- sum(adj[i, all.clusters[[4]]])
  vec5[i] <- sum(adj[i, all.clusters[[5]]])
  vec6[i] <- sum(adj[i, all.clusters[[6]]])
}
cbind(vec1, vec2, vec3, vec4, vec5, vec6) - previous.matrices$nxK.edge.counts$edge.counts.to

vec1 <- vec2 <- vec3 <- vec4 <- vec5 <- vec6 <- rep(0, 20);  
for(i in 1:20)
{
  vec1[i] <- sum(adj[all.clusters[[1]], i])
  vec2[i] <- sum(adj[all.clusters[[2]], i])
  vec3[i] <- sum(adj[all.clusters[[3]], i])
  vec4[i] <- sum(adj[all.clusters[[4]], i])
  vec5[i] <- sum(adj[all.clusters[[5]], i])
  vec6[i] <- sum(adj[all.clusters[[6]], i])
}
cbind(vec1, vec2, vec3, vec4, vec5, vec6) - previous.matrices$nxK.edge.counts$edge.counts.from
