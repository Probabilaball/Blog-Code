#Read in Data

d <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/Pitching Stabilization Points through Time/Historical Pitcher Data.csv',header=T)

#Define Functions for ML Estimation via the beta-binomial distribution


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

#startYear = the year to start the moving average
#endYear = the year to end the moving average
#nYears = number of years to include in moving average. 
#cutoff = minimum number of PA/AB required to be included in sample
#For Some statistics (SFs, for example), there is no data for early years. 
#Years start at 1900.


startYear = 2005
endYear = 2014
nYears = 6
cutoff = 100

#Create a vector of zeros to store stabiization points (the Ms)
#and the league mean talent level of each statistic (the mus)

MBHRFB = rep(0, length(startYear:endYear))

muHRFB = rep(0, length(startYear:endYear))

vHRFB = rep(0, length(startYear:endYear))
lHRFB = rep(0, length(startYear:endYear))
uHRFB = rep(0, length(startYear:endYear))

#Get vector of seasons to get the sample I want

season = d[,1]

#Begin Calculations

for(year in startYear:endYear) {

 #Only look at the number of years I'm interested in

 samp = (season <= year) & (season >= year - nYears +1)


 #HR/FB M Calculation

 x <- (d$HR)[samp]
 n <- (d$FB)[samp]

 s = (!is.na(x)) & (!is.na(n)) & (n >= cutoff) 

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MHRFB[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muHRFB[year-startYear+1] = ml$par[1]

 vPhi = v <- (-1/ml$par[2]^2)^2*ml$v[2,2]

 lHRFB[year-startYear+1] = (1-ml$par[2])/ml$par[2] - 1.96*sqrt(vPhi)
 uHRFB[year-startYear+1] = (1-ml$par[2])/ml$par[2] + 1.96*sqrt(vPhi)

 

}

#Make plots - use ggplot2 for time plots

library(ggplot2)


#Plots of stabilization Points over time (uses ggplot2)

qplot(startYear:endYear, MHRFB, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving BABIP Stabilization Point") + geom_line(col=2) + geom_line(aes(x = startYear:endYear, y = lHRFB), lty=2) +  geom_line(aes(x = startYear:endYear, y = uHRFB), lty=2)

#Plots of means over time (uses ggplot2)

qplot(startYear:endYear, muHRFB, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving BABIP Mean") + geom_line(col=2) 


