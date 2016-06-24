#Read in Data

d <- read.csv('C:/Users/Robert/Desktop/Multinomial-Dirichlet Empirical Bayes/2016 Hitting Data.csv', header=T)

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

sum(w*post.trout/sum(post.trout))



