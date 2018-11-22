

previous.matrices$KxK.edge.counts
previous.matrices$KxK.max.counts
previous.matrices$nxK.edge.counts$edge.counts.to
previous.matrices$nxK.edge.counts$edge.counts.from

## KxK edge count checks
for(i in 1:4)
{
  print(sum(adj[all.clusters[[i]], all.clusters[[i]]]))
}
sum(adj[all.clusters[[1]], all.clusters[[2]]])
sum(adj[all.clusters[[2]], all.clusters[[1]]])

sum(adj[all.clusters[[1]], all.clusters[[3]]])
sum(adj[all.clusters[[3]], all.clusters[[1]]])

sum(adj[all.clusters[[3]], all.clusters[[4]]])
sum(adj[all.clusters[[4]], all.clusters[[3]]])

sum(adj[all.clusters[[2]], all.clusters[[3]]])
sum(adj[all.clusters[[3]], all.clusters[[2]]])

## KxK max count checks
all.clusters
sapply(all.clusters, function(x) {MaxNodesWithinCluster(length(x), directed = TRUE)})
length(all.clusters[[3]]) * length(all.clusters[[2]])

## n x K matrices checks
vec1 <- rep(0, 20); vec2 <- rep(0, 20); vec3 <- rep(0, 20); vec4 <- rep(0, 20);
for(i in 1:20)
{
  vec1[i] <- sum(adj[i, all.clusters[[1]]])
  vec2[i] <- sum(adj[i, all.clusters[[2]]])
  vec3[i] <- sum(adj[i, all.clusters[[3]]])
  vec4[i] <- sum(adj[i, all.clusters[[4]]])
}
cbind(vec1, vec2, vec3, vec4) - previous.matrices$nxK.edge.counts$edge.counts.to

vec1 <- rep(0, 20); vec2 <- rep(0, 20); vec3 <- rep(0, 20);  vec4 <- rep(0, 20);
for(i in 1:20)
{
  vec1[i] <- sum(adj[all.clusters[[1]], i])
  vec2[i] <- sum(adj[all.clusters[[2]], i])
  vec3[i] <- sum(adj[all.clusters[[3]], i])
  vec4[i] <- sum(adj[all.clusters[[4]], i])
}
cbind(vec1, vec2, vec3, vec4) - previous.matrices$nxK.edge.counts$edge.counts.from