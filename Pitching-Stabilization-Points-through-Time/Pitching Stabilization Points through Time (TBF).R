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


startYear = 1925
endYear = 2014
nYears = 6
cutoff = 400

#Create a vector of zeros to store stabiization points (the Ms)
#and the league mean talent level of each statistic (the mus)

MSO = rep(0, length(startYear:endYear))
MOBP = rep(0, length(startYear:endYear))
MH = rep(0, length(startYear:endYear))
MOBB = rep(0, length(startYear:endYear))
MHBP = rep(0, length(startYear:endYear))
MHR = rep(0, length(startYear:endYear))


muSO = rep(0, length(startYear:endYear))
muOBP = rep(0, length(startYear:endYear))
muH = rep(0, length(startYear:endYear))
muBB = rep(0, length(startYear:endYear))
muHBP = rep(0, length(startYear:endYear))
muHR = rep(0, length(startYear:endYear))

#Get vector of seasons to get the sample I want

season = d[,1]

#Begin Calculations

for(year in startYear:endYear) {

 #Only look at the number of years I'm interested in

 samp = (season <= year) & (season >= year - nYears +1)


 #SO M Calculation

 x <- (d$SO)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MSO[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muSO[year-startYear+1] = ml$par[1]

 #OBP M Calculation

 x <- (d$H + d$BB + d$HBP)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MOBP[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muOBP[year-startYear+1] = ml$par[1]

 #BB M Calculation

 x <- (d$BB)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MBB[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muBB[year-startYear+1] = ml$par[1]

 #HBP M Calculation

 x <- (d$HBP)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MHBP[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muHBP[year-startYear+1] = ml$par[1]

 #H M Calculation

 x <- (d$H)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MH[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muH[year-startYear+1] = ml$par[1]

 #HR M Calculation

 x <- (d$HR)[samp]
 n <- (d$TBF)[samp]

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart = sum(x)/sum(n)
 N = length(x)
 s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
 MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
 phiStart = 1/(MStart+1)

 ml = mlBetaBinom(c(muStart,phiStart),x,n)
 
 MHR[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 muHR[year-startYear+1] = ml$par[1]


}

#Make plots - use ggplot2 for time plots

library(ggplot2)

#Look at all stabilization points for all years at once - neat picture, but a jumbled mess


#Plots of stabilization Points over time (uses ggplot2)

qplot(startYear:endYear, MSO, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving SO Rate Stabilization Point") + geom_line(col=2)
qplot(startYear:endYear, MOBP, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving OBP Rate Stabilization Point") + geom_line(col=2)
qplot(startYear:endYear, MH, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving H Rate Stabilization Point") + geom_line(col=2)
qplot(startYear:endYear, MBB, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving BB Rate Stabilization Point") + geom_line(col=2)
qplot(startYear:endYear, MHBP, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving HBP Rate Stabilization Point") + geom_line(col=2)
qplot(startYear:endYear, MHR, xlab = 'Year', ylab = 'Estimated Stabilization Point', main = "Six Year Moving HR Rate Stabilization Point") + geom_line(col=2)


#Plots of means over time (uses ggplot2)

qplot(startYear:endYear, muSO, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving SO Rate Mean") + geom_line(col=2)
qplot(startYear:endYear, muOBP, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving OBP Rate Mean") + geom_line(col=2)
qplot(startYear:endYear, muH, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving H Rate Mean") + geom_line(col=2)
qplot(startYear:endYear, muBB, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving BB Rate Mean") + geom_line(col=2)
qplot(startYear:endYear, muHBP, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving HBP Rate Mean") + geom_line(col=2)
qplot(startYear:endYear, muHR, xlab = 'Year', ylab = 'Estimated Mean', main = "Six Year Moving HR Rate Mean") + geom_line(col=2)
