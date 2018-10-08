
## Results
last.iter <- 27201
plot(ts(num.clusters[1:last.iter]), xlab = "iters", ylab = "number of clusters")
IAT(num.clusters)


no.as.iat <- c(44.23607, 43.09678, 44.5097, 38.78507)
as.iat <- c(43.74782, 38.73524, 44.15076, 36.74659)
num.particles <- c(5, 10, 20, 50)


plot(num.particles, no.as.iat, type = 'b', ylab = "Integrated Autocorrelation Time", 
     xlab = "Number of particles")
lines(num.particles, as.iat, type = 'b', lty = 3)
legend(38, 45, c("Without A.S.", "With A.S."), lty = c(1, 3), cex = 0.8)

