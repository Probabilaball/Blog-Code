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


#Define functions to create Covariance matrices


cov.matrix <- function(pars) {
  
   k <- length(pars)
   m <- matrix(rep(0,k*k),k,k)

   a0 <- sum(pars)

   for(i in 1:k) {
    for(j in 1:k) {
      
      if(i == j) m[i,j] = pars[i]*(a0-pars[i])/(a0^2*(a0 + 1))
      if(i != j) m[i,j] = -pars[i]*pars[j]/(a0^2*(a0 + 1))
    }
   }

   return(m)

}

cov.pred.matrix <- function(pars,n.new) {
  
   k <- length(pars)
   m <- matrix(rep(0,k*k),k,k)

   a0 <- sum(pars)

   for(i in 1:k) {
    for(j in 1:k) {
      
      if(i == j) m[i,j] = 1/n.new*pars[i]/a0*(1-pars[i]/a0)*(n.new+a0)/(1 + a0)
      if(i != j) m[i,j] = -1/n.new*pars[i]*pars[j]/a0^2*(n.new+a0)/(1+a0)
    }
   }

   return(m)

}



#Specific example: Mike Trout in 2013



#Define array of wOBA weights

w <- c(.89,1.27,1.62,2.10,0.69,0.72,0)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,100,9,407)

#Calculate dirichlet posterior for Mike Trout 2013
#And take a weighted average of expectations

post.trout <- x.trout + alpha

woba.trout <- sum(w*post.trout/sum(post.trout))

#Estimate posterior predictive wOBA distribution by Normal Approximation

v <- t(w)%*%cov.matrix(post.trout)%*%w

c(woba.trout - 1.96*sqrt(v), woba.trout + 1.96*sqrt(v))



#Estimate predictive wOBA distribution by normal approximation

v <- t(w)%*%cov.pred.matrix(post.trout, sum(x.trout))%*%w

c(woba.trout - 1.96*sqrt(v), woba.trout + 1.96*sqrt(v))

