#Read in Data

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/Regressed-Confidence-Intervals-for-wOBA/2016 Hitting Data.csv', header=T)

#Select only players with at least 300 events
#Event is defined as AB + BB - IBB + SF + HBP

n <- d$AB + d$BB - d$IBB + d$SF + d$HBP 

samp <- (n >= 300)

d <- d[samp,]
n <- n[samp]

#Create matrix of counts of outcomes
#Outcomes are single, double, triple, home run, walk, and hit by pitch

x <- as.matrix(cbind(d$X1B, d$X2B, d$X3B, d$HR, d$BB, d$HBP))

#Append final column consisting of all other outcomes to event

x <- cbind(x, n - apply(x,1,sum))

#Fit the multinomial-dirichlet model using the dirmult package

library(dirmult)

dirmult.fit <- dirmult(x)

alpha <- dirmult.fit$gamma


#Specific example: Mike Trout in 2013



#Define array of wOBA weights

w <- c(.89,1.27,1.62,2.10,0.69,0.72,0)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,110,9,397)

#Calculate dirichlet posterior for Mike Trout 2013
#And take a weighted average of expectations

post.trout <- x.trout + alpha



#Estimate posterior wOBA distribution by simulation

n.iter <- 500000

sim.data <- rdirichlet(n.iter, post.trout)
sim.woba <- w[1]*sim.data[,1] + w[2]*sim.data[,2] + w[3]*sim.data[,3] + w[4]*sim.data[,4] + w[5]*sim.data[,5] + w[6]*sim.data[,6]

hist(sim.woba, freq = F, main = 'Histogram of simulated posterior wOBA values', xlab = 'Simulated wOBA')

quantile(sim.woba, c(.025,.975))


#Check normality

qqnorm(sim.woba)
qqline(sim.woba)




#Estimate posterior predictive wOBA distribution by simulation

sim.pred <- matrix(rep(0, n.iter * 7), n.iter, 7)

for(i in 1:n.iter) sim.pred[i,] <- rmultinom(1, sum(x.trout) , sim.data[i,])

sim.pred.woba <- rep(0, n.iter)

for(i in 1:n.iter) sim.pred.woba[i] <- sum(w*sim.pred[i,]/sum(sim.pred[i,]))

hist(sim.pred.woba, freq = F, main = 'Histogram of simulated posterior predictive wOBA values', xlab = 'Simulated wOBA')

quantile(sim.pred.woba, c(.025,.975))


#Check normality

qqnorm(sim.pred.woba)
qqline(sim.pred.woba)

