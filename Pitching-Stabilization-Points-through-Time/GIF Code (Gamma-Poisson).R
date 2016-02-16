#Read in Data

d <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/Pitching Stabilization Points through Time/Historical Pitcher Data.csv',header=T)

#Define Functions for ML Estimation via the beta-binomial distribution


#Define Functions for ML Estimation via the beta-binomial distribution


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

#startYear = the year to start the moving average
#endYear = the year to end the moving average
#nYears = number of years to include in moving average. 
#cutoff = minimum number of PA/AB required to be included in sample
#For Some statistics (SFs, for example),
#there is no data for early years. Years start at 1920.


startYear = 1920
endYear = 2014
nYears = 6
cutoff = 80

#Create a vector of zeros to store stabiization points (the Ms)
#and the league mean talent level of each statistic (the mus)

M = rep(0, length(startYear:endYear))
mu = rep(0, length(startYear:endYear))

#Create a vector of lower and upper bounds for when I make the plots
#I run the code once to find the lower and upper x and y limits
#Then replace those in the code that plots all the images. 

xMin = 1
xMax = 0
yMax = 0

#Get vector of seasons to look at

season = d[,1]

#Begin Calculations

for(year in startYear:endYear) {

 #Only look at the number of years I'm interested in

 samp = (season <= year) & (season >= year - nYears +1)

 n <- d$IP[samp]
 n <- round(n) + (n-round(n))/.3
 x <- (d$TBF[samp]-3*n)
 n <- (d$IP)[samp]
 n <- round(n) + (n-round(n))/.3

 s = (!is.na(x)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 muStart <- mean(x/n)
 MStart <- mean((x/n)^2)/var(x/n)
 phiStart <- 1/(1+MStart)

 ml = mlGammaPoisson(c(muStart,phiStart),x,n)
 
 M[year-startYear+1] = (1-ml$par[2])/ml$par[2]
 mu[year-startYear+1] = ml$par[1]

 #This is me keeping track of what the smallest and largest
 #x and y values that both the observed data
 #and the estimated talent distribution take

 xMin = min(xMin, x/n)
 xMax = max(xMax, x/n)
 yMax = max(yMax, max(hist(x/n, plot=F)$density))
 yMax = max(yMax, max(dgamma((ml$par[1]*(1-ml$par[2])/ml$par[2]-1)/((1-ml$par[2])/ml$par[2]), ml$par[1]*(1-ml$par[2])/ml$par[2], (1-ml$par[2])/ml$par[2])))

 #This is the code to actually make the plots for the gif
 #Remove commenting when you figure out the appropriate values
 #for the x and y limits and number of breaks

 plotName <- paste(year, ".jpg")
 jpeg(plotName)
 hist(x/n, freq=F, xlab = "XBF Rate", ylab = "Density", main = year, xlim = c(0.65,2.15), ylim = c(0, 3.8), breaks = seq(0.65,2.15, by = 0.05)) 
 curve(dgamma(x, ml$par[1]*(1-ml$par[2])/ml$par[2], (1-ml$par[2])/ml$par[2]), add=T, lty=2)
 dev.off()

 #The plots get saved in my documents folder. Maybe different on yours.
 #I use a freeware gif program to make a gif.

}






