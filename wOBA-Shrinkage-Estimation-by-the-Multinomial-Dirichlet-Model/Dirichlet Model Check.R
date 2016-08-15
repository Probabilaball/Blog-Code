#Read in Data

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/wOBA-Shrinkage-Estimation-by-the-Multinomial-Dirichlet-Model/2016 Hitting Data.csv', header=T)

#Select only players with at least 300 events
#Event is defined as AB + BB - IBB + SF + HBP

n <- d$AB + d$BB - d$IBB + d$SF + d$HBP 

samp <- (n >= 300)

d <- d[samp,]
n <- n[samp]

#Create matrix of counts of outcomes
#Outcomes are single, double, triple, home run, walk, and hit by pitch

x <- as.matrix(cbind(d$X1B, d$X2B, d$X3B, d$HR, d$BB - d$IBB, d$HBP))

#Append final column consisting of all other outcomes to event

x <- cbind(x, n - apply(x,1,sum))


#Fit the multinomial-dirichlet model using the dirmult package

library(dirmult)

dirmult.fit <- dirmult(x)

alpha <- dirmult.fit$gamma

#Simulate thetas (sets of talent levels) from fitted Dirichlet distribution

theta <- rdirichlet(length(n), alpha)

#Simulate observed counts of outcomes using n's from the data set

sim.x <- matrix(rep(0, length(n) * 7), length(n), 7)

for(i in 1:length(n)) sim.x[i,] <- rmultinom(1, n[i] , theta[i,])


#Compare means and standard deviations of observed and simulated data

data.frame('Empirical' = round(apply(x/n, 2, mean),4), 'Simulated' = round(apply(sim.x/n, 2, mean),4))

data.frame('Empirical' = round(apply(x/n, 2, sd),4), 'Simulated' = round(apply(sim.x/n, 2, sd),4))

#Side-by-side histograms of real and simulated data
#for each component

par(mfrow = c(2,7))

for(i in 1:7) hist(x[,i]/n, freq=F)
for(i in 1:7) hist(sim.x[,i]/n, freq=F)