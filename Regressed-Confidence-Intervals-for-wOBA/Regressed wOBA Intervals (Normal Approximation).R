#Define functions to create Covariance matrices
#for standard and prediction intervals
#given a set of Dirichlet parameters pars
#and potentially a number of new plate appearances n.new

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

#Define array of Dirichlet prior parameters

alpha <- c(34.30, 10.44, 1.16, 5.74, 16.29, 1.96, 144.51)

#Define array of counts of events for Mike Trout 2013

x.trout <- c(115,39,9,27,100,9,407)

#Calculate parameters of Dirichlet posterior
#distribution for Mike Trout 2013

post.trout <- x.trout + alpha

#Calculate wOBA for Mike Trout 2013

woba.trout <- sum(w*post.trout/sum(post.trout))

#Estimate posterior predictive wOBA distribution by Normal Approximation

v <- t(w)%*%cov.matrix(post.trout)%*%w

c(woba.trout - 1.96*sqrt(v), woba.trout + 1.96*sqrt(v))



#Estimate predictive wOBA distribution by normal approximation

v <- t(w)%*%cov.pred.matrix(post.trout, sum(x.trout))%*%w

c(woba.trout - 1.96*sqrt(v), woba.trout + 1.96*sqrt(v))

