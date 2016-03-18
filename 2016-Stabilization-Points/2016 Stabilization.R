#Read in Data

h <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/2016-Stabilization-Points/2016 Hitting Data.csv',header=T)
p <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/2016-Stabilization-Points/2016 Pitching Data.csv',header=T)

#h = hitter data
#p = pitcher data

#Define functions


invHess <- function(m) {

 d <- m[1,1]*m[2,2] - m[1,2]*m[2,1]
 m2 <- matrix(c(m[2,2], -m[1,2], -m[2,1], m[1,1]), 2,2,byrow=T)
 return(m2/d)

}


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

plot.stabilization.p <- function(M, se, lower.bound = 0.25, upper.bound = 0.75) {

   #Plots estimated n as a function of stabilization level

   p <- seq(lower.bound, upper.bound, by = .001)

   n.est <- p/(1-p)*M

   l <- n.est - 1.96*p/(1-p)*se
   u <- n.est + 1.96*p/(1-p)*se

   plot(p, n.est, type = 'l', xlab = "Stabilization Level", ylab = "Required Number of Events", main = "Estimated Number of Events Required\nfor Given Stabilization Level")
   points(p,l, type= 'l', lty=2)
   points(p,u, type = 'l', lty=2)

}

plot.stabilization.n <- function(M, se, lower.n = 1, upper.n = 600) {

   #Plots estimated stabilization  level as a function of n

   n <- seq(lower.n, upper.n,by = 1)

   p.est <- n/(n + M)

   l <- p.est - 1.96*n/(n+M)^2*se
   u <- p.est + 1.96*n/(n+M)^2*se

   plot(n,p.est, type = 'l', xlab = "Number of Events", ylab = "Estimated Stabilization Level", main = "Estimated Stabilization Level\nfor Given Number of Events",ylim = c(0,max(u)))
   points(n,l, type= 'l', lty=2)
   points(n,u, type = 'l', lty=2)

}

stabilize <- function(x, n, cutoff, model) {

 x <- x[n >= cutoff]
 n <- n[n >= cutoff]

 #Beta-binomial model

 if(model == 1) {

   muStart = sum(x)/sum(n)
   N = length(x)
   s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
   MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
   phiStart = 1/(MStart+1)

 
   ml = mlBetaBinom(c(muStart,phiStart),x,n)
   mu = ml$par[1]
   phi = ml$par[2]
   M = (1-phi)/phi
  
   se <- sqrt((-1/phi^2)^2*ml$v[2,2])

   plot.stabilization.n(M, se,upper.n=max(n))

   return(list('Parameters' = c(mu, M), "Standard.Error" = se, "Model" = "Beta-Binomial"))

  }
  
  if(model == 2) {

   muStart <- mean(x/n)
   MStart <- mean((x/n)^2)/var(x/n)
   phiStart <- 1/(1+MStart)


   ml = mlGammaPoisson(c(muStart,phiStart),x,n)
   mu = ml$par[1]
   phi = ml$par[2]
   M = (1-phi)/phi

   se <- sqrt((-1/phi^2)^2*ml$v[2,2])

   plot.stabilization.n(M, se,upper.n=max(n))

   return(list('Parameters' = c(mu, M), "Standard.Error" = se, "Model" = "Gamma-Poisson"))

  }

}
 
  
#Tell R to show lots of digits.
#Otherwise you end up with lots of 
#output in scientific notation.

options('scipen' = 999)

#OFFENSIVE STATS

#On-Base Percentage

stabilize(h$H + h$BB + h$HBP, h$PA, cutoff = 300, 1)   

#Batting Average

stabilize(h$H,h$AB, cutoff = 300, 1)

#Strikeout Rate

stabilize(h$SO, h$PA, cutoff=300,1)

#Walk Rate

stabilize(h$BB - h$IBB, h$PA - h$IBB, cutoff=300,1)

#Singles Rate

stabilize(h$X1B,h$PA, cutoff = 300, 1)

#Doubles Rate

stabilize(h$X2B,h$PA, cutoff = 300, 1)

#Triples Rate

stabilize(h$X3B,h$PA, cutoff = 300, 1)

#XBH Rate

stabilize(h$X2B + h$X3B,h$PA, cutoff = 300, 1)

#HR Rate

stabilize(h$HR, h$PA, cutoff=300,1)

#HBP Rate

stabilize(h$HBP, h$PA, cutoff=300,1)

#BABIP

stabilize(h$H - h$HR, h$AB - h$SO - h$HR + h$SF, cutoff = 300, 1)



#PITCHING STATS

#Correct innings pitched

IP <- p$IP
IP <- round(IP) + (IP-round(IP))/.3

#BABIP

stabilize(p$H - p$HR, p$GB + p$FB + p$LD, cutoff = 300,1)

#GB Rate

stabilize(p$GB, p$GB + p$FB + p$LD, cutoff = 300,1)

#FB Rate

stabilize(p$FB, p$GB + p$FB + p$LD, cutoff = 300,1)

#LD Rate

stabilize(p$LD, p$GB + p$FB + p$LD, cutoff = 300,1)

#HR/FB Rate

stabilize(p$HR,p$FB, cutoff = 100, 1)

#SO Rate

stabilize(p$SO, p$TBF, cutoff = 400, 1)

#HR Rate

stabilize(p$HR, p$TBF, cutoff = 400, 1)

#BB Rate

stabilize(p$BB - p$IBB, p$TBF - p$IBB, cutoff = 400, 1)

#HBP Rate

stabilize(p$HBP, p$TBF, cutoff = 400, 1)

#H Rate

stabilize(p$H, p$TBF, cutoff = 100, 1)

#OBP

stabilize(p$H + p$BB + p$HBP, p$TBF, cutoff=400,1)

#WHIP, ER Rate, and XBF use Gamma-Poisson because the stat is not
#a binary outcome which takes only 0 or 1 per event (happens or does
#not happen), it is a count which can take 0, 1, 2, 3, etc.

#Whip

stabilize(p$BB + p$H, IP, cutoff = 80,2)

#ER Rate

stabilize(p$ER, IP, cutoff = 80,2)

#Extra batters faced

stabilize(p$TBF - 3*IP, IP, cutoff = 80, 2)

