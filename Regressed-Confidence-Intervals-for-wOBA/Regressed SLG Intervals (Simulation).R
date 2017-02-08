#Load the dirmult package for simulation

library(dirmult)

#Specific example: Mike Trout in 2013

#Define array of SLG weights

w <- c(1,2,3,4,0)

#Define array of Dirichlet prior parameters

alpha <- c(42.44, 12.86, 1.38, 7.07, 176.12)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,399)

#Calculate dirichlet posterior for Mike Trout 2013
#And take a weighted average of expectations

post.trout <- x.trout + alpha



#Estimate posterior SLG distribution by simulation

#set the seed for consistent results

set.seed(5)


#Take n.iter simulations

n.iter <- 500000

#Simulate true talent levels for the posterior for Mike Trout in 2013
#and calculate wOBA for each

sim.data <- rdirichlet(n.iter, post.trout)
sim.slg <- w[1]*sim.data[,1] + w[2]*sim.data[,2] + w[3]*sim.data[,3] + w[4]*sim.data[,4]

#Take quantiles for a 95% credible interval

slg.ci <- quantile(sim.slg, c(.025,.975))

slg.ci

#Show histogram of values

hist(sim.slg, freq = F,  main = 'Histogram of simulated posterior SLG values\nfor Mike  Trout in 2013', xlab = 'Simulated SLG')

abline(v = slg.ci, lty = 2)

#Check normality using a quantile-quantile plot

qqnorm(sim.slg)
qqline(sim.slg)




#Estimate posterior predictive SLG distribution by simulation

set.seed(5)

sim.pred <- matrix(rep(0, n.iter * 5), n.iter, 5)

for(i in 1:n.iter) sim.pred[i,] <- rmultinom(1, sum(x.trout) , sim.data[i,])

sim.pred.slg <- rep(0, n.iter)

for(i in 1:n.iter) sim.pred.slg[i] <- sum(w*sim.pred[i,]/sum(sim.pred[i,]))

pred.slg.ci <- quantile(sim.pred.slg, c(.025,.975))


hist(sim.pred.slg, freq = F, main = 'Histogram of simulated posterior predictive wOBA \nvalues for Mike Trout in 2013', xlab = 'Simulated wOBA')
abline(v = pred.slg.ci, lty = 2)

#Check normality using quantile-quantile plot

qqnorm(sim.pred.slg)
qqline(sim.pred.slg)

