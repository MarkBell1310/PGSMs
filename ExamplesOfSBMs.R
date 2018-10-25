
#****************************************************
#
#****** Examples of stochastic block models  ********
#
#****************************************************

## Generating SBMs

sbm <- sample_sbm(n = 20,
                  pref.matrix = matrix(c(0.99, 0, 0, 0.1), c(2, 2)),
                  block.sizes = c(10, 10), directed = FALSE, loops = FALSE); plot(sbm)
num.clusters <- 4
sbm <- sample_sbm(n = 20,
                  pref.matrix = forceSymmetric(matrix(runif(num.clusters^2),
                                                      (c(num.clusters, num.clusters)))),
                  block.sizes = c(6, 4, 7, 3), directed = FALSE, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20,
                  pref.matrix = forceSymmetric(
                    matrix(c(1, 0.01, 0.01, 1), c(2, 2))),
                  block.sizes = c(10, 10), directed = FALSE, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20, pref.matrix = diag(5),
                  block.sizes = rep(4, 5), directed = FALSE, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20, pref.matrix = diag(4),
                  block.sizes = rep(5, 4), directed = FALSE, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = 20, pref.matrix = 0.99 * diag(1),
                  block.sizes = 20, directed = FALSE, loops = FALSE); plot(sbm)

## Starting clusters

start.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
                       c(2, 6, 8, 10, 11, 17))
start.clusters <- list(1:10, 11:num.nodes)
start.clusters <- list(1:num.nodes)


sbm <- sample_sbm(n = 50, pref.matrix = diag(5),
                  block.sizes = rep(10, 5), directed = FALSE, loops = FALSE); plot(sbm)

# generate SBM
num.nodes <- 20 # no. nodes
K <- 4  # no. clusters
#pref.matrix <- diag(K) # Bernoulli rates (K x K matrix)
#pref.matrix = forceSymmetric(matrix(runif(K^2),(c(K, K))))
pref.matrix = forceSymmetric(matrix(rbeta(K^2, 2, 2),(c(K, K))))
block.sizes <- rep(num.nodes/K, K) # no. nodes in each cluster (K length vector)
directed <- FALSE
sbm <- sample_sbm(num.nodes, pref.matrix, block.sizes, directed = directed, loops = FALSE); plot(sbm)
start.clusters <- list(1:num.nodes) # list of clusters
adj <- as_adj(sbm)

# ERGM stats
y <- network(as.matrix(adj), directed=FALSE)
summary(y ~ edges + kstar(2)); plot(y)

# generate SBM
num.nodes <- 20 # no. nodes
directed <- FALSE
start.clusters <- list(1:num.nodes) # list of clusters
#start.clusters <- list(1:25, 26:50)
#start.clusters <- list(1:4, 5:8, 9:12, 13:16, 17:20)
#start.clusters <- list(1:10, 11:20, 21:30, 31:40, 41:50)
#start.clusters <- start.clusters <- list(1:5, 6:10, 11:15, 16:20, 21:25, 26:30, 31:35, 36:40, 41:45, 46:50)
#start.clusters <- sapply(1:num.nodes, function(x){list(x)})
# start.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
#                        c(2, 6, 8, 10, 11, 17))
#sbm <- sample_sbm(n = num.nodes, pref.matrix = diag(5),
#                 block.sizes = rep(4, 5), directed = directed, loops = FALSE); plot(sbm)
sbm <- sample_sbm(n = num.nodes, pref.matrix = diag(5),
                  block.sizes = rep(4, 5), directed = FALSE, loops = FALSE); plot(sbm)
adj <- as_adj(sbm)