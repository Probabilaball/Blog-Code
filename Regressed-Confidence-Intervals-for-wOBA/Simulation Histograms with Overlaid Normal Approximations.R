#For standard interval

hist(sim.woba, freq = F, main = 'Histogram of simulated posterior wOBA values\nfor Mike Trout in 2013', xlab = 'Simulated wOBA')
curve(dnorm(x, woba.trout, sqrt(v)), add=T)

#For prediction interval

hist(sim.pred.woba, freq = F, main = 'Histogram of simulated posterior predictive wOBA\nvalues for Mike Trout in 2013', xlab = 'Simulated wOBA')
curve(dnorm(x, woba.trout, sqrt(v)), add=T)
