#Read in data

d <- read.csv('C:/[Your file path]/obp data.csv',header=T)

#Functions for estimating M

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
 dmuphi = dphimu = sum(1/phi^2*(-digamma(-mu/phi+mu+1/phi-1)+digamma(mu/phi-mu)-(mu-1)*(phi-1)/phi*trigamma(-mu/phi+mu+1/phi-1)+mu*(1/phi-1) *trigamma(mu/phi-mu)+digamma(n-x+mu-mu/phi+1/phi-1)+(mu-1)*(phi-1)/phi*trigamma(n-x+mu-mu/phi+1/phi-1)-digamma(x-mu+mu/phi) +mu*(phi-1)/phi*trigamma(x-mu+mu/phi)))
 dphiphi = sum(1/phi^4*(-mu^2*trigamma(mu/phi-mu)+(mu-1)^2*(-trigamma(-mu/phi+mu+1/phi-1))+2*(mu-1)*phi*digamma(-mu/phi+mu+1/phi-1)-2*mu*phi*digamma(mu/phi-mu)+(mu-1)^2*trigamma(n-x+mu-mu/phi+1/phi-1)-2*(mu-1)*phi*digamma(n-x+mu-mu/phi+1/phi-1)-2*phi*digamma (n+1/phi-1)-trigamma(n+1/phi-1)+mu^2*trigamma(x-mu+mu/phi)+2*mu*phi*digamma(x-mu+mu/phi)+2*phi*digamma(1/phi-1)+trigamma(1/phi  -1)))

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

  h <- -invHess(hess(x,n, par))
  return(list(par = par, 'var' = h))

}

#Lower = Smallest number of PA/AB/Whatever to consider
#Upper = Largest number of PA/AB/Whatever to consider

lower = 150
upper = 600

#Vectors of zeroes to store things

M_ = rep(0, upper-lower+1)
mu_ = rep(0, upper-lower+1)

l_ = rep(0, upper-lower+1)
u_ = rep(0, upper-lower+1)

#Start loop for calculations

for(cutoff in lower:upper) {
  
 #Statistic of interest - change this to change what 
 #Statistic you are looking at. Currently set to OBP.

 x <- d$H + d$BB + d$HPB
 n <- d$PA

 x <- x[n >= cutoff]
 n <- n[n >= cutoff]

 #Starting values for ML algorithm
 
 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 #Maximize and store the results

 ml = mlBetaBinom(c(muStart,phiStart),x,n) 

 mu_[cutoff-lower+1] = ml$par[1]
 M_[cutoff-lower+1] = (1-ml$par[2])/ml$par[2]


 #Calculate and store lower and upper confidence bounds

 vPhi <- ml$v[2,2]
 v <- (-1/ml$par[2]^2)^2*vPhi

 l_[cutoff-lower+1] = (1-ml$par[2])/ml$par[2]-1.96*sqrt(v)
 u_[cutoff-lower+1] = (1-ml$par[2])/ml$par[2]+1.96*sqrt(v)


}

#Plot estimated M as a function of cutoff point, with confidence bounds

plot(lower:upper,M_, xlab = 'Cutoff Point', ylab = 'Estimated M', type = 'l', ylim = c(min(l_)-1, max(u_)+1))
points(lower:upper, l_, type = 'l', lty=2)
points(lower:upper, u_, type = 'l', lty = 2)
