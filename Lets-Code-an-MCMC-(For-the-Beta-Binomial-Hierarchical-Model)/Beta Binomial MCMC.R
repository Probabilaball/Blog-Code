#Read in data

d <- read.csv('C:/Users/Robert/Desktop/Lets Code an MCMC (For the Beta-Binomial Hierarchical Model)/2015 Batting Average Data.csv')

x <- d$H[d$AB >= 300]
n <- d$AB[d$AB >= 300]

hist(x/n, main = "Histogram of Batting Averages\nin 2015 MLB Season (Min 300 AB)", xlab = "Observed Batting Average")

#Set the seed so I get the same results each time I run the code.

set.seed(5)

#MCMC Function


betaBin.mcmc <- function(x, n, mu.start, phi.start, burn.in = 1000, n.draws = 5000, sigma.mu = 0.005, sigma.phi = 0.001) {

 #Input variables: x, n, mu.start, phi.start burn.in, n.iter, sigmaMu, sigmaPhi
 #x is the number of successes in n trials. These are the data being input.
 #x and n have no default value. mu.start and phi.start are starting values
 #for the chain values of mu and phi - different starting values should be 
 #chosen for several different chains to check for convergence.
 #burn.in is the burn-in period of the mcmc chain, defaulted to 1000. 
 #n.draws is the number of draws from the mcmc chain after burn-in, defaulted 
 #to 5000. sigma.mu and sigma.phi are standard deviations
 #for the normal candidate distributions in the Metropolis-Hastings
 #step for mu and phi. These are defaulted to 0.005 and 0.001, but 
 #will have to be changed for different data sets. The goal is 
 #to manipulate these so that the acceptance rates for mu and phi are
 #each at around 40%.


 #cond = conditional likelihood function for beta-binomial
 #input parameter values phi, mu and data x,n
 #Return sum of log-likelihood (l) and log prior (p)

 cond = function(phi, mu, x, n) {
  	 l = sum(lbeta(mu*(1-phi)/phi + x, (1-mu)*(1-phi)/phi+n-x) - lbeta(mu*(1-phi)/phi, (1-mu)*(1-phi)/phi))
       p = dbeta(mu,0.5,0.5,log=T)+dbeta(phi,0.5,0.5,log=T)
       return(l + p)
 }


 #Create vectors/matrices of zeroes to store MCMC draws

 phi = rep(0, burn.in + n.draws)
 mu = rep(0, burn.in + n.draws)
 theta = matrix(rep(0, length(n)*(burn.in + n.draws)), length(n), (burn.in + n.draws))

 #Create variables to store acceptance values

 acceptance.mu = 0
 acceptance.phi = 0

 #Define starting values for mu and phi
 #Using mu.start and phi.start entered
 #as values into the function

 mu[1] = mu.start
 phi[1] = phi.start

 #Simulate starting values for thetas using initial mu and phi values

 for(j in 1:length(x)) {
    theta[j, 1] = rbeta(1, mu[1]*(1-phi)[1]/phi[1] + x[j], (1-phi)[1]/phi[1]*(1-mu[1]) + n[j] - x[j])
 }

 #Begin MCMC chain. This chain will run for the length of the burn-in period
 #specified plus the number of iterations requested

 for(i in 2:(burn.in + n.draws)) {

   #Set current chain values equal to previous chain values

   phi[i] = phi[i-1]
   mu[i] = mu[i-1]

   #Metropolis-Hastings step for mu

   #Draw a candidate observation from a normal distribution
   #centered at the previous observation in the chain
   #with a pre-specified variance.

   cand = rnorm(1,mu[i-1],sigma.mu)

   #Check if candidate is between 0 and 1. If not, discard it.
  
   if((cand > 0) & (cand < 1)) {

	  

  	  lpo = cond(phi[i-1],mu[i-1],x, n)
 	  lpn = cond(phi[i-1],cand,x, n)

	  #Draw an observation from a uniform(0,1)

 	  u = runif(1)


        if((lpn - lpo) > log(u)) {
		mu[i] = cand
		acceptance.mu = acceptance.mu+1
	  }
  }

  #Metropolis-Hastings step for phi

  #Draw a candidate observation from a normal distribution
  #centered at the previous observation in the chain
  #with a pre-specified variance.

  cand = rnorm(1,phi[i-1],sigma.phi)

  #Check if candidate is between 0 and 1. If not, discard it.
  
  if( (cand > 0) & (cand < 1)) {
         
 	lpo = cond(phi[i-1],mu[i-1],x, n)
	lpn = cond(cand,mu[i-1],x, n)


	#Draw an observation from a uniform(0,1)

 	u = runif(1)

      if((lpn - lpo) > log(u)) {
		phi[i] = cand
		acceptance.phi = acceptance.phi + 1
	} 


  }

  #Gibbs sampling step for theta
  #For each player, simulate a batting average from a beta distribution
  #With parameters alpha = x + mu*(1-phi)/phi and beta = n - x + (1-mu)*(1-phi)/phi
  #where mu and phi are the current values in the chain.

  for(j in 1:length(n)) {
   theta[j, i] = rbeta(1, (1-phi[i])/phi[i]*mu[i] + x[j], (1-phi[i])/phi[i]*(1-mu[i]) + n[j] - x[j])
  } 

 } #End MCMC Loop

 #Throw out the burn-in iterations

 mu <- mu[(burn.in + 1):(burn.in + n.draws)]
 phi <- phi[(burn.in + 1):(burn.in + n.draws)]
 theta <- theta[,(burn.in + 1):(burn.in + n.draws)]

 #Return chain for each parameter and acceptance rates for MH step parameters

 return(list(mu = mu, phi = phi, theta = theta, acceptance = c(acceptance.mu/(burn.in + n.draws), acceptance.phi/(burn.in + n.draws))))

}

#Use multiple chains to and check that all chains are converging to the same
#stationary distribution after the burn-in. Use multiple starting values to
#help check.

chain.1 <- betaBin.mcmc(x,n, 0.265, 0.002)
chain.2 <- betaBin.mcmc(x,n, 0.5, 0.1)
chain.3 <- betaBin.mcmc(x,n, 0.100, 0.0001)

#Visually inspect the chains for mu, phi, and one theta to see if
#all chains are converging to the same stationary distribution. A big blob of
#colors is good. One or more chains separate from the rest is bad.

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")
matplot(data.frame(chain.1$phi, chain.2$phi, chain.3$phi), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Phi", main = "MCMC Chains for Phi")
matplot(data.frame(chain.1$theta[1,], chain.2$theta[1,], chain.3$theta[1,]), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Theta 1", main = "MCMC Chains for Theta 1")

#Check acceptance rates for mu and phi. Around 40% is good.

chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

#Inference

#Combine all three chains into one for inference

mu <- c(chain.1$mu, chain.2$mu, chain.3$mu)
phi <- c(chain.1$phi, chain.2$phi, chain.3$phi)
theta <- cbind(chain.1$theta, chain.2$theta, chain.3$theta)


#Mean, standard deviation, and 95% intervals for mu

hist(mu, xlab = 'Mu', freq=F)
mean(mu)
sd(mu)
quantile(mu,c(.025,.975))

#Mean, standard deviation, and 95% intervals for phi

hist(phi, xlab = 'Phi', freq=F)
mean(phi)
sd(phi)
quantile(phi,c(.025,.975))

#Mean, standard deviation, and 95% intervals for theta, the
#batting average for each player in the sample

names <- d$Name[d$AB >= 300]

hist(theta[1,], xlab = 'Batting Averge', main = 'Posterior Distribution of Batting\nAverage for Brye Harper',freq=F)

data.frame("Player" = names, "Mean" = apply(theta,1,mean), "SD" = apply(theta,1,sd), "Lower 95" = apply(theta,1,quantile, 0.025), "Upper 95" = apply(theta,1,quantile, 0.975))

