#Read in data

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/Lets-Code-an-MCMC-(For-the-Beta-Binomial-Hierarchical-Model)/2015 Batting Average Data.csv')

#Cut sample to only players with 300 or more AB in the season 

x <- d$H[d$AB >= 300]
n <- d$AB[d$AB >= 300]


#This is a method of moments estimator for the parameters of the beta-binomial
#distribution or two-stage model. This estimator is originally described in
#Kleinman, J. (1973), "Proportions with Extraneous Variance: Single and Independent Sample,"
#Journal of the American Statistical Association, 68(341), 46 - 54.

mm.est <- function(x, n) {

 #Input: x, vector of successes, and n, vector of trials
 
 #Initialize starting estimates of true binomial means
 #As raw proportions

 p.i <- x/n

 #Initialize weights as sample sizes
 #And create a number to store it

 w.i <- n

 #diff = difference, keep track of how sum of squared difference
 #between new and old weights. When difference gets lower than a certain
 #tolerance, stop.
 #iter = count of iterations. If iterations gets too high, stop and throw an error

 diff <- 10000
 iter <- 0

 while(diff > 10^(-6)) { 

  iter <- iter + 1

  #w is the sum of weights

  w <- sum(w.i)

  #p is the estimated mean of the underlying beta distribution

  p <- sum(w.i*p.i)/w

  #s is a weighted estimate of the second moment

  s <- sum(w.i*(p.i - p)^2)

  #phi is the estimated dispersion parameter of the underlying beta distribution

  phi <- (s-p*(1-p)*sum(w.i/n *(1-w.i/w)))/(p*(1-p)*(sum(w.i*(1-w.i/w))-sum(w.i/n*(1-w.i/w))))

  #Re-calculate weights and squared sum difference between new and old weights

  w.i.old <- w.i

  w.i <-  n/(1+phi*(n-1))

  diff <- sum((w.i.old - w.i)^2)

  if(iter > 1000) stop("Algorithm failed to converge")

  }

  #Convert answers to the traditional alpha, beta format

  return(list('alpha' = p*(1-phi)/phi, 'beta' = (1-p)*(1-phi)/phi))

}

#Estimate the parameters of the underlying beta distribution

betabin.par <- mm.est(x,n)

alpha <- betabin.par$alpha
beta <- betabin.par$beta

#Mean and 95% interval for Bryce Harper's batting average


bryce.alpha <- x[1] + alpha
bryce.beta <- n[1] - x[1] + beta

bryce.alpha/(bryce.alpha + bryce.beta)

qbeta(c(.025,.975), bryce.alpha, bryce.beta)

curve(dbeta(x, bryce.alpha, bryce.beta), xlim = c(0.24, 0.36), xlab = 'Batting Averge', main = 'Posterior Distribution of Batting\nAverage for Bryce Harper')


#Mean, standard deviation, and 95% intervals for theta, the
#batting average for each player in the sample

names <- d$Name[d$AB >= 300]

data.frame("Player" = names, "Mean" = (x + alpha)/(n + alpha + beta), "Lower 95" = qbeta(0.025, x + alpha, n - x + beta), "Upper 95" = qbeta(0.975, x + alpha, n - x + beta))




