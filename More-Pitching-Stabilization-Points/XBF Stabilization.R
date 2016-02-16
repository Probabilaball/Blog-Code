#Read in Data

data <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/More Pitcher Stabilization Points/Pitcher Data.csv',header=T)

#Define functions


#This is my code for finding the maximum likelihood estimates via 
#the Newton-Rhapson Algorithm. I've defined the log likelihood function
#even though I don't use it. Grad is the matrix of first derivatives of the
#log likelihood function and hess is the matrix of second derivatives of the
#log likelihood function.

gradBB = function(x,n,par) {

  mu = par[1]
  phi = par[2]

  dmu = sum(1/phi*(phi-1)*(digamma(mu*(1/phi-1))-digamma((mu-1)*(phi-1)/phi)+digamma(n-x+mu-mu/phi+1/phi-1)-digamma(x+mu*(1/phi-1))))
  dphi = sum((1-mu)/phi^2*digamma(-mu/phi+mu+1/phi-1)+mu/phi^2*digamma(mu/phi-mu)+(mu-1)/phi^2*digamma(n-x+mu-mu/phi+1/phi-1)+1/phi^2*digamma(n+1/phi-1)-mu/phi^2*digamma(x-mu+mu/phi)-1/phi^2*digamma(1/phi-1))

  return(matrix(c(dmu,dphi),2,1))

}

hessBB = function(x,n,par) {

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

   diff <- -invHess(hessBB(x,n, par))%*%gradBB(x,n, par)
   par <- par + diff

  }

  v <- -invHess(hessBB(x,n, par))
  return(list(par = par, 'var' = v))

}

gradGP = function(x,n,par) {

  mu = par[1]
  phi = par[2]

  dmu = sum((phi-1)/phi*(digamma(mu*(1/phi-1)) + log((n-1)*phi+1) - digamma(x + mu*(1/phi-1)) - log(1-phi)))
  dphi = sum((mu*(1-phi)+x*phi)/(phi^2*((n-1)*phi+1))+mu/phi^2*(log(((n-1)*phi+1)/phi) + digamma(-mu*(phi-1)/phi) - digamma((mu*(1-phi)+x*phi)/phi) - log(-(phi-1)/phi)-1))
  return(matrix(c(dmu,dphi),2,1))

}


hessGP = function(x,n,par) {

 mu = par[1]
 phi = par[2]

 dmumu = sum((1-phi)^2/phi^2*trigamma(x + mu*(1/phi-1)) - (1-phi)^2/phi^2*trigamma(mu*(1/phi-1)))
 dphiphi = sum(-mu^2/phi^4*trigamma(mu/phi-mu) + 3*mu/phi^3  + 2*mu/phi^3*log(1-phi) + 2*mu/phi^3*log(1/phi) - 2*mu/phi^3*digamma(mu/phi-mu) + mu/((1-phi)*phi^2) + mu/(phi^5*(n + (1-phi)/phi)^2) - 4*mu/(phi^4 *(n + (1-phi)/phi)) - mu/(phi^4*(n + (1-phi)/phi)^2) + 2*mu/(phi^3*(n + (1-phi)/phi)) - 2*mu/phi^3*log(n + (1-phi)/phi) + x/(phi^4*(n + (1-phi)/phi)^2) - 2*x/(phi^3*(n + (1-phi)/phi)) + mu^2/phi^4*trigamma(x-mu + mu/phi) + 2*mu/phi^3*digamma(x-mu + mu/phi))
 dmuphi = dphimu = sum(mu/phi^3*trigamma(mu*(1-phi)/phi) + 1/phi^2*digamma(mu*(1-phi)/phi) - mu/phi^2*trigamma(mu*(1-phi)/phi) + 1/(phi^3*(n + (1-phi)/phi)) - 1/(phi^2*(n + (1-phi)/phi)) + log(n + (1-phi)/phi)/phi^2 - mu/phi^3*trigamma(x + mu*(1-phi)/phi) - 1/phi^2 *digamma(x + mu*(1-phi)/phi) + mu/phi^2*trigamma(x + mu*(1-phi)/phi)-1/phi^2 - log((1-phi)/phi)/phi^2)

 return(matrix(c(dmumu,dmuphi,dphimu,dphiphi),2,2))

}

mlGammaPoisson <- function(par, x, n) {

  diff = 100

  while( sum(diff^2) > 10^(-10)) {

   diff <- -invHess(hessGP(x,n, par))%*%gradGP(x,n, par)
   par <- par + diff

   par <- abs(par)
   if((par[2] <= 0) | (par[2] >= 1)) par[2] = runif(1)

  }

  h <- -invHess(hessGP(x,n, par))
  return(list(par = par, 'var' = h))

}

n <- data$IP
n <- round(n) + (n-round(n))/.3
x <- data$TBF - 3*n

n <- round(n) + (n-round(n))/.3

x <- x[n >= 80]
n <- n[n >= 80]

#muStart = sum(x)/sum(n)
#N = length(x)
#s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
#MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
#phiStart = 1/(MStart+1)


muStart <- mean(x/n)
MStart <- mean((x/n)^2)/var(x/n)
phiStart <- 1/(1+MStart)


#Maximize it

ml = mlGammaPoisson(c(muStart,phiStart),x,n)
mu = ml$par[1]
phi = ml$par[2]
M = (1-phi)/phi

#yMax <- dbeta((mu*M-1)/(M-2), mu*M, (1-mu)*M)
yMax <- dgamma((mu*M-1)/M, mu*M, M)

hist(x/n,freq=F, ylim = c(0, yMax), xlab = "XBF/IP",ylab = "Density", main = "Observed XBF Rate with True Talent Distribution")
#curve(dbeta(x, mu*M, (1-mu)*M), add=T, lty = 2)
curve(dgamma(x, mu*M, M), add=T, lty=2)

v <- (-1/phi^2)^2*ml$v[2,2]

M

sqrt(v)

M-1.96*sqrt(v)
M+1.96*sqrt(v)


#CI Plot - define vector of stabilization points

p <- seq(0.4, 0.8, by = .001)

#Calculate central CI values

c <- p/(1-p)*M

#Calculate lower and upper 95% confidence bounds

l <- c - 1.96*p/(1-p)*sqrt(v)
u <- c + 1.96*p/(1-p)*sqrt(v)

#Plot the estimated values with confidence bounds

plot(p,c, type = 'l', xlab = "Stabilization Level", ylab = "Innings Pitched Required", main = "Innings Pitched Required for XBF Stabilization, 2009-2014")
points(p, l, type= 'l', lty=2)
points(p, u, type = 'l', lty=2)



