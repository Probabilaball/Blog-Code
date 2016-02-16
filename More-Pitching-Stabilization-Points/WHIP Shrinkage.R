#Read in the data

d <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/WHIP Stabilization by the Gamma-Poisson Model/Pitcher Data.csv',header=T)

#x is the sum number of hits + walks per pitcher
#n is the number of innings pitched

x <- d$H
n <- d$IP

#Do a little bit of extra manipulation to get the IP correct
#because motherfuckers still write 120.33 IP as 120.1 for some dumb reason

n <- round(n) + (n-round(n))/.3

#Only take pitchers with at least 80 IP

x <- x[n >= 80]
n <- n[n >= 80]


#Functions for ML Estimation

#Gradient = returns matrix of first derivatives of log-likelihood

grad = function(x,n,par) {

  mu = par[1]
  phi = par[2]

  dmu = sum((phi-1)/phi*(digamma(mu*(1/phi-1)) + log((n-1)*phi+1) - digamma(x + mu*(1/phi-1)) - log(1-phi)))
  dphi = sum((mu*(1-phi)+x*phi)/(phi^2*((n-1)*phi+1))+mu/phi^2*(log(((n-1)*phi+1)/phi) + digamma(-mu*(phi-1)/phi) - digamma((mu*(1-phi)+x*phi)/phi) - log(-(phi-1)/phi)-1))
  return(matrix(c(dmu,dphi),2,1))

}

#Hessian = returns matrix of second derivatives of log-likelihood function

hess = function(x,n,par) {

 mu = par[1]
 phi = par[2]

 dmumu = sum((1-phi)^2/phi^2*trigamma(x + mu*(1/phi-1)) - (1-phi)^2/phi^2*trigamma(mu*(1/phi-1)))
 dphiphi = sum(-mu^2/phi^4*trigamma(mu/phi-mu) + 3*mu/phi^3  + 2*mu/phi^3*log(1-phi) + 2*mu/phi^3*log(1/phi) - 2*mu/phi^3*digamma(mu/phi-mu) + mu/((1-phi)*phi^2) + mu/(phi^5*(n + (1-phi)/phi)^2) - 4*mu/(phi^4 *(n + (1-phi)/phi)) - mu/(phi^4*(n + (1-phi)/phi)^2) + 2*mu/(phi^3*(n + (1-phi)/phi)) - 2*mu/phi^3*log(n + (1-phi)/phi) + x/(phi^4*(n + (1-phi)/phi)^2) - 2*x/(phi^3*(n + (1-phi)/phi)) + mu^2/phi^4*trigamma(x-mu + mu/phi) + 2*mu/phi^3*digamma(x-mu + mu/phi))
 dmuphi = dphimu = sum(mu/phi^3*trigamma(mu*(1-phi)/phi) + 1/phi^2*digamma(mu*(1-phi)/phi) - mu/phi^2*trigamma(mu*(1-phi)/phi) + 1/(phi^3*(n + (1-phi)/phi)) - 1/(phi^2*(n + (1-phi)/phi)) + log(n + (1-phi)/phi)/phi^2 - mu/phi^3*trigamma(x + mu*(1-phi)/phi) - 1/phi^2 *digamma(x + mu*(1-phi)/phi) + mu/phi^2*trigamma(x + mu*(1-phi)/phi)-1/phi^2 - log((1-phi)/phi)/phi^2)

 return(matrix(c(dmumu,dmuphi,dphimu,dphiphi),2,2))

}

#invHess = returns the inverse of a 2x2 matrix m

invHess <- function(m) {

 d <- m[1,1]*m[2,2] - m[1,2]*m[2,1]
 m2 <- matrix(c(m[2,2], -m[1,2], -m[2,1], m[1,1]), 2,2,byrow=T)
 return(m2/d)

}

#mlGammaPoisson = performs ML estimation by Newton-Raphson algorithm
#Input is starting points par and data x and n. Output is list containing
#ml estimates and variance/covariance matrix

mlGammaPoisson <- function(par, x, n) {

  diff = 100

  while( sum(diff^2) > 10^(-10)) {

   diff <- -invHess(hess(x,n, par))%*%grad(x,n, par)
   par <- par + diff

   #I was having a problem where sometimes the updated value of the parameters
   #would produce a negative phi estimate, which would crash the program.
   #I solved this by taking the absolute value of the parameters after 
   #each step (basically restarting it at a new spot if the values are 
   #negative). I'm sure there is a better way to solve this problem.

   par <- abs(par)
   if((par[2] <= 0) | (par[2] >= 1)) par[2] = runif(1)

  }

  h <- -invHess(hess(x,n, par))
  return(list(par = par, 'var' = h))

}

#Calculate moment estimators as starting values

muStart <- mean(x/n)
kStart <- mean((x/n)^2)/var(x/n)
phiStart <- 1/(1+kStart)

#Get ML Estimates and store them in proper variables

ml = mlGammaPoisson(c(muStart, phiStart),x,n)

mu = ml$par[1]
phi = ml$par[2]
k = (1-ml$par[2])/ml$par[2]

#Get variance of K by delta method

vPhi <- ml$v[2,2]
v <- (-1/phi^2)^2*vPhi

#Plot observed WHIP values and overlay estimated true talent distribution

hist(x/n, freq=F, ylim = c(0,3), xlab = "WHIP", ylab = "Density", main = "Observed WHIP with True Talent Distribution")
curve(dgamma(x, mu*k, k), add=T, lty = 2)


#CI Plot - define vector of stabilization points

p <- seq(0.4, 0.8, by = .001)

#Calculate central CI values

c <- p/(1-p)*k

#Calculate lower and upper 95% confidence bounds

l <- c - 1.96*p/(1-p)*sqrt(v)
u <- c + 1.96*p/(1-p)*sqrt(v)

#Plot the estimated values with confidence bounds

plot(p,c, type = 'l', xlab = "Stabilization Level", ylab = "Innings Pitched Required", main = "IP Required for WHIP Stabilization, 2009-2014, min 80 IP")
points(p, l, type= 'l', lty=2)
points(p, u, type = 'l', lty=2)



