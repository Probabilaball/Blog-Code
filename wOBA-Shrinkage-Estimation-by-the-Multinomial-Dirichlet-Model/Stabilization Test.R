#Set number of simulated wOBA values

n.sims <- 50000

#Read in hitting data, 2010 - 2015

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/wOBA-Shrinkage-Estimation-by-the-Multinomial-Dirichlet-Model/2016 Hitting Data.csv', header=T)

#Fit Model for players with at least 300 PA

samp <- (d$PA  >= 300)

d <- d[samp,]


n <- d$PA
x <- as.matrix(cbind(d$X1B, d$X2B, d$X3B, d$HR, d$BB - d$IBB, d$HBP))
x <- cbind(x, n - apply(x,1,sum))



library(dirmult)
res <- dirmult(cbind(x, n - apply(x,1,sum)))

#Store dirichlet parameters in a vector

par <- res$gamma

#n.PA = alpha0 = sum of Dirichlet paramaters = stabilization point

n.PA <- round(sum(par))

#Store wOBA weights in a vector

w <- c(.89,1.27,1.62,2.10,0.69,0.72,0)

#Simulate thetas 
#These are the sets of joint talent levels for the categories

theta <- rdirichlet(n.sims, par)

#For each talent level, simulate two sets of counts of events
#for each set of thetas, both based on a sample of size n.pa

x.1 <- matrix(rep(0, n.sims*length(par)), n.sims, length(par))
for(i in 1:n.sims) x.1[i,] <- t(rmultinom(1, n.PA, theta[i,]))

x.2 <- matrix(rep(0, n.sims*length(par)), n.sims, length(par))
for(i in 1:n.sims) x.2[i,] <- t(rmultinom(1, n.PA, theta[i,]))

#Calculate the wOBA for each of the simulated counts of events

woba.1 <- rep(0, n.sims)
woba.2 <- rep(0, n.sims)

for(i in 1:n.sims) woba.1[i] <- sum(w*x.1[i,]/n.PA)
for(i in 1:n.sims) woba.2[i] <- sum(w*x.2[i,]/n.PA)

#Plot the two simulated sets of wOBA and calculate the correlation

plot(woba.1, woba.2)
cor(woba.1, woba.2)

