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

#This is a function for using the bootstrap with the method of moments estimator above
#to correct the length of empirical Bayesian intervals in the beta-binomial model
#More details on the parametric bootstrap correction method for empirical Bayesian intervals
#Can be found in Laird, N., and Louis, T. (1987), "Empirical Bayes Confidence Intervals Based on 
#Bootstrap Samples," Journal of the American Statistical Association, 82(399), 739-750.


mm.bootstrap <- function(x, n, nBootstraps = 5000, confidence = 0.95) {

 #Inputs are x, a vector of counts, in n observations total
 #pre-specificed inputs are the number of bootstraps to perform
 #nBootstraps and the confidence level

 #Get estimated parameters of the underlying distribution alpha and beta

 par <- mm.est(x, n) 

 alpha <- par$alpha
 beta <- par$beta

 #Construct naive (uncorrected) empirical Bayesian intervals
 #alphaEB and betaEB are parameters of the posterior distributions for each observation
 #naiveL and naiveU are uncorrected lower and upper central bounds of the posterior. Default
 #confidence level is set at 95%

 alphaEB <- x + alpha
 betaEB <- n - x + beta
  
 naiveL <- qbeta((1-confidence)/2, alphaEB, betaEB)
 naiveU <- qbeta(confidence + (1-confidence)/2, alphaEB, betaEB)

 #alphaBoot and betaBoot are vectors to store the parametric
 #bootstrap estimates

 alphaBoot <- rep(0, nBootstraps)
 betaBoot <- rep(0, nBootstraps)

 for(i in 1:nBootstraps) {

  #Generate data from two-stage model
  #Using parameter estimates from method of moments estimator

  #An extra check is added to ensure a vector of all x/n = 0
  #or all x/n = 1 is not generated

  xBoot <- rep(0, length(x)) 
  
  while( sum(xBoot) == 0 | sum(xBoot) == sum(n)) {

   theta <- rbeta(length(x), alpha, beta)
   xBoot <- rbinom(length(x), n, theta)

  }

  #Store bootstrapped parameter estimates in appropriate vectors
  #Run a try() in case a data set is generated
  #that causes the algorithm to fail to converge. These should
  #be rare but something to watch out for. In the event this happens
  #the bootstrapped parameters from the previous iteration will be used.

  parBoot <- try(mm.est(xBoot, n),TRUE)

  #A character string should be thrown back if an error occurs

  if(!is.character(parBoot[1])) {

   alphaBoot[i] <- parBoot$alpha
   betaBoot[i] <- parBoot$beta

  } else{

   alphaBoot[i] <- alphaBoot[i-1]
   betaBoot[i] <- betaBoot[i-1]

  }

 }

 #Define function to give average quantile value over a set of bootstrapped
 #parameters, minus a goal. 

 p <- function(x, quantile, alphaBoot, betaBoot) mean(pbeta(x, alphaBoot, betaBoot),na.rm=T) - quantile
  
 #bootL and bootU are upper and lower interval values for the bootstrap-corrected
 #empirical Bayesian intervals

 bootL <- rep(0, length(x))
 bootU <- rep(0, length(x))

 #Solve for bootstrap intervals

 for(i in 1:length(x)) {

  alphaBootEB <- x[i] + alphaBoot
  betaBootEB <- n[i] - x[i] + betaBoot

  bootL[i] <- uniroot(p, interval=c(0,1), quantile = (1-confidence)/2, alphaBootEB, betaBootEB)$root
  bootU[i] <- uniroot(p, interval=c(0,1), quantile = confidence + (1-confidence)/2, alphaBootEB, betaBootEB)$root
 
 }

 #Return naive intervals and bootstrap intervals

 return(list('Naive' = data.frame(naiveL, naiveU), 'Bootstrap' = data.frame(bootL, bootU)))

}


   
#One word of caution - it is common (even likely) that the bootstrapped parameters
#alpha and beta will contain some values. These are ignored when calculating the
#quantile, averaging over the bootstrapped values. This means that even though
#you have 10000 bootstraps, you may be throwing some away and averaging over
#fewer. These tend to be extreme value and will not have much of an effect
#on the intervals. This is also where the warnings are coming from - calculating
#pbeta(x, alphaBootEB, betaBootEB) when sometimes alphaBootEB or betaBootEB is negative

#Example 1: Data used by Morris, Efron in empirical Bayes example
#Originally found in 
#Efron, B and Morris, C (1975), "Data Analysis Using Stein's Estimator and its Generalizations",
#Journal of the American Statistical Association, 70(350), 311 - 319.
#x is the number of hits in n at-bats

x <- c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)
x <- round(x*45)
n <- rep(45, 18)


mm.est(x,n)
mm.bootstrap(x,n)

#Example 2: Simulated data from a beta distribution
#with mean = 0.1 and relatively large variance
#Unequal sample sizes used 

n <- rpois(10, 10)
theta <- rbeta(length(n), .1*4, (1-.1)*4)
x <- rbinom(length(n), n, theta)

mm.est(x,n)
mm.bootstrap(x,n)
