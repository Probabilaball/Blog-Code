#Read in Data

data <- read.csv('C:/[Your file path here]/obp data.csv',header=T)

#Define Random variables. x is the observed on-base events in n plate appearances
#Currently set to OBP, but this may be changed here (plot titles will still show OBP, however)

x <- data$H + data$BB + data$HBP
n <- data$PA

#Only select those with greater than 300 plate appearances

x <- x[n >= 300]
n <- n[n >= 300]

#This is my code for finding the maximum likelihood estimates via 
#the Newton-Rhapson Algorithm. I've defined the log likelihood function
#even though I don't use it. Grad is the matrix of first derivatives of the
#log likelihood function and hess is the matrix of second derivatives of the
#log likelihood function.

logLike = function(param, x, n) {
 mu = param[1]
 phi =  param[2]
 N = length(x) 
 t = sum(lbeta((1-phi)/phi*mu + x, (1-phi)/phi*(1-mu)+n-x)) - N*lbeta((1-phi)/phi*mu, (1-phi)/phi*(1-mu)) 
 return(t)
}

grad = function(x,n,par) {

  mu = par[1]
  phi = par[2]

  dmu = sum(1/phi*(phi-1)*(digamma(mu*(1/phi-1))-digamma((mu-1)*(phi-1)/phi)+digamma(n-x+mu-mu/phi+1/phi-1)-digamma(x+mu*(1/phi-1))))
  dphi = sum((1-mu)/phi^2*digamma(-mu/phi+mu+1/phi-1)+mu/phi^2*digamma(mu/phi-mu)+(mu-1)/phi^2*digamma(n-x+mu-mu/phi+1/phi-1)+1/phi^2*digamma(n+1/phi-1)-mu/phi^2*digamma(x-mu+mu/phi)-1/phi^2*digamma(1/phi-1))

  return(matrix(c(dmu,dphi),2,1))

}

hess = function(x,n,par) {

 mu = par[1]
 phi = par[2]

 dmumu = sum(-(phi-1)^2/phi^2*trigamma(mu*(1/phi-1))-(phi-1)^2/phi^2*trigamma((mu-1)*(phi-1)/phi)+(phi-1)^2/phi^2*trigamma(n-x+mu- mu/phi+1/phi-1)+(phi-1)^2/phi^2*trigamma(x+mu*(1/phi-1)))
 dphiphi = sum(1/phi^4*(-mu^2*trigamma(mu/phi-mu)+(mu-1)^2*(-trigamma(-mu/phi+mu+1/phi-1))+2*(mu-1)*phi*digamma(-mu/phi+mu+1/phi-1)-2*mu*phi*digamma(mu/phi-mu)+(mu-1)^2*trigamma(n-x+mu-mu/phi+1/phi-1)-2*(mu-1)*phi*digamma(n-x+mu-mu/phi+1/phi-1)-2*phi*digamma (n+1/phi-1)-trigamma(n+1/phi-1)+mu^2*trigamma(x-mu+mu/phi)+2*mu*phi*digamma(x-mu+mu/phi)+2*phi*digamma(1/phi-1)+trigamma(1/phi-1)))
 dmuphi = dphimu = sum(1/phi^2*(-digamma(-mu/phi+mu+1/phi-1)+digamma(mu/phi-mu)-(mu-1)*(phi-1)/phi*trigamma(-mu/phi+mu+1/phi-1)+mu*(1/phi-1) *trigamma(mu/phi-mu)+digamma(n-x+mu-mu/phi+1/phi-1)+(mu-1)*(phi-1)/phi*trigamma(n-x+mu-mu/phi+1/phi-1)-digamma(x-mu+mu/phi) +mu*(phi-1)/phi*trigamma(x-mu+mu/phi)))

 return(matrix(c(dmumu,dmuphi,dphimu,dphiphi),2,2))

}

invHess <- function(m) {

 d <- m[1,1]*m[2,2] - m[1,2]*m[2,1]
 m2 <- matrix(c(m[2,2], -m[1,2], -m[2,1], m[1,1]), 2,2,byrow=T)
 return(m2/d)

}

mlBetaBinom <- function(par, x, n) {

  diff = 100

  while( sum(diff^2) > 10^(-10)) {

   diff <- -invHess(hess(x,n, par))%*%grad(x,n, par)
   par <- par + diff

  }

  v <- -invHess(hess(x,n, par))
  return(list(par = par, 'var' = v))

}


#Get starting values for maximization algorithm from moment estimators

muStart = sum(x)/sum(n)
N = length(x)
s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
phiStart = 1/(MStart+1)

#Maximize it

ml = mlBetaBinom(c(muStart,phiStart),x,n)
mu = ml$par[1]
phi = ml$par[2]
M = (1-phi)/phi

#Make a histogram of observed OBP values and overlay the estimated talent distribution 

hist(x/n,freq=F, xlab = "OBP", ylim = c(0,14),ylab = "Probability", main = "Distribution of Observed OBP with True Talent Distribution")
curve(dbeta(x, mu*M, (1-mu)*M), add=T, lty = 2)


#CI Plot - define vector of stabilization points

p <- seq(0.4, 0.8, by = .001)

#Calculate central CI values

c <- p/(1-p)*M

#vPhi is the estimate of the variance of Phi = 1/(M+1)

vPhi <- ml$v[2,2]

#Transform variance to the scale of M 
#And Calculate 95% confidence interval for each p 


v <- (-1/phi^2)^2*vPhi

l <- c - 1.96*p/(1-p)*sqrt(v)
u <- c + 1.96*p/(1-p)*sqrt(v)

#Plot the estimated values with confidence bounds

plot(p,c, type = 'l', xlab = "Stabilization Level", ylab = "Plate Appearances Required", main = "PA Required for OBP Stabilization, 2009-2014")
points(p, l, type= 'l', lty=2)
points(p, u, type = 'l', lty=2)
