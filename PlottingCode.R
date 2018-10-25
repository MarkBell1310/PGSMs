
## Bayes Afternoon meeting results

plot(num.particles, no.as.iat, type = 'b', ylab = "Integrated Autocorrelation Time", 
     xlab = "Number of particles")
lines(num.particles, as.iat, type = 'b', lty = 3)
legend(38, 45, c("Without A.S.", "With A.S."), lty = c(1, 3), cex = 0.8)

algorithm <- c("Gibbs", "PGSM (20)", "AS (20)", "PGSM & Gibbs")
correct.K <- c(8.9, 21.2, 17.6, 23.1)
IAT <- c(79.7, 43.4, 39.5, 11.1)

#correct.K <- c(8.9, 21.2, 22.7, 17.6, 18.1, 23.1)
#IAT <- c(79.7, 43.4, 41.2, 39.5, 41.0, 11.1)


# Correct K
plot(correct.K, xaxt = "n", xlab='Algorithm', type = 'b', ylab = "% correct K")
axis(1, at=1:4, labels = c("Gibbs", "PGSM (20)", "AS (20)", "PGSM & Gibbs"))

# IAT
plot(IAT, xaxt = "n", xlab='Algorithm', type = 'b', ylab = "IAT")
axis(1, at=1:4, labels = c("Gibbs", "PGSM (20)", "AS (20)", "PGSM & Gibbs"))













