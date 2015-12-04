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


#Example 1: Data used by Morris, Efron in empirical Bayes example
#Originally found in 
#Efron, B and Morris, C (1975), "Data Analysis Using Stein's Estimator and its Generalizations",
#Journal of the American Statistical Association, 70(350), 311 - 319.
#x is the number of hits in n at-bats

x <- c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)
x <- round(x*45)
n <- rep(45, 18)


mm.est(x,n)

#Example 2: Simulated data from a beta distribution
#Unequal sample sizes used 

n <- rpois(10, 50)
theta <- rbeta(10, 2, 2)
x <- rbinom(10, n, theta)

mm.est(x,n)


