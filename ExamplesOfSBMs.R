
#****************************************************
#
#****** Examples of stochastic block models  ********
#
#****************************************************

## Generating SBMs

sbm <- sample_sbm(n = 20,
                  pref.matrix = matrix(c(0.99, 0, 0, 0.1), c(2, 2)),
                  block.sizes = c(10, 10), directed = FALSE, loops = FALSE); plot(sbm)
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
sbm <- sample_sbm(n = 20, pref.matrix = 0.99 * diag(1),
                  block.sizes = 20, directed = FALSE, loops = FALSE); plot(sbm)

## Starting clusters

start.clusters <- list(c(18, 14, 3, 5), c(12, 16, 20), c(1, 4, 19), c(7, 9, 13, 15),
                       c(2, 6, 8, 10, 11, 17))
start.clusters <- list(1:10, 11:num.nodes)
start.clusters <- list(1:num.nodes)
start.clusters <- list(1:3, 4:7, 8:10, 11:13, 14:17, 18:20) 