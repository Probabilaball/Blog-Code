#Read in Data

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/Regressed-Confidence-Intervals-for-wOBA/2016 Hitting Data.csv', header=T)

#Select only players with at least 300 events
#Event is defined as AB

n <- d$AB  

samp <- (n >= 300)

d <- d[samp,]
n <- n[samp]

#Create matrix of counts of outcomes
#Outcomes are single, double, triple, home run

x <- as.matrix(cbind(d$X1B, d$X2B, d$X3B, d$HR))

#Append final column consisting of all other outcomes to event

x <- cbind(x, n - apply(x,1,sum))

#Fit the multinomial-dirichlet model using the dirmult package

library(dirmult)

dirmult.fit <- dirmult(x)

alpha <- dirmult.fit$gamma


#Specific example: Mike Trout in 2013



#Define array of SLG weights

w <- c(1,2,3,4,0)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,399)

#Calculate dirichlet posterior for Mike Trout 2013
#And take a weighted average of expectations

post.trout <- x.trout + alpha



#Estimate posterior SLG distribution by simulation

n.iter <- 500000

sim.data <- rdirichlet(n.iter, post.trout)
sim.slg <- w[1]*sim.data[,1] + w[2]*sim.data[,2] + w[3]*sim.data[,3] + w[4]*sim.data[,4]

hist(sim.slg, freq = F, main = 'Histogram of simulated SLG values', xlab = 'Simulated SLG')

quantile(sim.slg, c(.025,.975))


#Check normality

qqnorm(sim.slg)
qqline(sim.slg)




#Estimate posterior predictive SLG distribution by simulation

sim.pred <- matrix(rep(0, n.iter * 5), n.iter, 5)

for(i in 1:n.iter) sim.pred[i,] <- rmultinom(1, sum(x.trout) , sim.data[i,])

sim.pred.slg <- rep(0, n.iter)

for(i in 1:n.iter) sim.pred.slg[i] <- sum(w*sim.pred[i,]/sum(sim.pred[i,]))

hist(sim.pred.slg, freq = F, main = 'Histogram of simulated predictive SLG values', xlab = 'Simulated SLG')

quantile(sim.pred.slg, c(.025,.975))


#Check normality

qqnorm(sim.pred.slg)
qqline(sim.pred.slg)

