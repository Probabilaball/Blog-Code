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

sum(w*post.trout/sum(post.trout))



