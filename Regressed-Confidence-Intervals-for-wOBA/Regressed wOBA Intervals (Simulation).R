#Load the dirmult package for simulation

library(dirmult)

#Specific example: Mike Trout in 2013

#Define array of wOBA weights

w <- c(.89,1.27,1.62,2.10,0.69,0.72,0)

#Define array of Dirichlet prior parameters

alpha <- c(34.30, 10.44, 1.16, 5.74, 16.29, 1.96, 144.51)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,100,9,407)

#Calculate parameters of Dirichlet posterior
#distribution for Mike Trout 2013

post.trout <- x.trout + alpha

#Estimate posterior wOBA distribution by simulation


#set the seed for consistent results

set.seed(5)

#Take n.iter simulations

n.iter <- 500000

#Simulate true talent levels for the posterior for Mike Trout in 2013
#and calculate wOBA for each

sim.data <- rdirichlet(n.iter, post.trout)
sim.woba <- w[1]*sim.data[,1] + w[2]*sim.data[,2] + w[3]*sim.data[,3] + w[4]*sim.data[,4] + w[5]*sim.data[,5] + w[6]*sim.data[,6]


#Take quantiles for credible interval

woba.ci <- quantile(sim.woba, c(.025,.975))

woba.ci

#Show histogram of values

hist(sim.woba, freq = F, main = 'Histogram of simulated posterior wOBA values\nfor Mike Trout in 2013', xlab = 'Simulated wOBA')
abline(v = woba.ci, lty = 2)

#Check normality using a quantile-quantile plot

qqnorm(sim.woba)
qqline(sim.woba)




#Estimate posterior predictive wOBA distribution by simulation

set.seed(5)

sim.pred <- matrix(rep(0, n.iter * 7), n.iter, 7)

for(i in 1:n.iter) sim.pred[i,] <- rmultinom(1, sum(x.trout) , sim.data[i,])

sim.pred.woba <- rep(0, n.iter)

for(i in 1:n.iter) sim.pred.woba[i] <- sum(w*sim.pred[i,]/sum(sim.pred[i,]))


pred.ci <- quantile(sim.pred.woba, c(.025,.975))

pred.ci

hist(sim.pred.woba, freq = F, main = 'Histogram of simulated posterior predictive wOBA\nvalues for Mike Trout in 2013', xlab = 'Simulated wOBA')
abline(v = pred.ci, lty = 2)

#Check normality using quantile-quantile plot

qqnorm(sim.pred.woba)
qqline(sim.pred.woba)

