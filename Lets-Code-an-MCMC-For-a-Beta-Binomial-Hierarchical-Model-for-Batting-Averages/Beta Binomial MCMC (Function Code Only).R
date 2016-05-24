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
 #each at roughly around 40%.


 #m = log of function for beta-binomial that is proportional to full posterior
 #distribution for mu and phi (it does not include normalizing constants).
 #Input constant parameter values mu, phi and data vectors x,n
 #Return sum of log-likelihood (l) and log prior (p)
 #this is the function I call m in my blog post.

 m  = function(mu, phi, x, n) {
       N = length(x)
       l = sum(lbeta(mu*(1-phi)/phi + x, (1-mu)*(1-phi)/phi+n-x)) - N*lbeta(mu*(1-phi)/phi, (1-mu)*(1-phi)/phi)
       p = -0.5*log(mu) - 0.5*log(1-mu) - 0.5*log(phi) - 0.5*log(1-phi)
       return(l + p)
 }


 #Create vectors/matrices of zeroes to store MCMC draws

 phi = rep(0, burn.in + n.draws)
 mu = rep(0, burn.in + n.draws)
 theta = matrix(rep(0, length(n)*(burn.in + n.draws)), length(n), (burn.in + n.draws))

 #Create counter variables to store frequency of acceptance in MH steps

 acceptance.mu = 0
 acceptance.phi = 0

 #Define starting values for mu and phi
 #Using mu.start and phi.start entered
 #as values into the function

 mu[1] = mu.start
 phi[1] = phi.start

 #Simulate starting values for thetas using initial mu and phi values

 for(i in 1:length(x)) {
    theta[i, 1] = rbeta(1, mu[1]*(1-phi)[1]/phi[1] + x[i], (1-phi)[1]/phi[1]*(1-mu[1]) + n[i] - x[i])
 }

 #Begin MCMC chain. This chain will run for the length of the burn-in period
 #specified plus the number of iterations requested

 for(j in 2:(burn.in + n.draws)) {

   #Set current chain values equal to previous chain values

   phi[j] = phi[j-1]
   mu[j] = mu[j-1]

   #Metropolis-Hastings step for mu

   #Draw a candidate observation from a normal distribution
   #centered at the previous observation in the chain
   #with a pre-specified variance.

   cand = rnorm(1, mu[j-1], sigma.mu)

   #Check if candidate is between 0 and 1. If not, discard it.
  
   if((cand > 0) & (cand < 1)) {

	  #m.old = Proportional log posterior function at all previous values
	  #m.new = Proportional log posterior function at all previous values
	  #	     except for candidate mu value.

  	  m.old = m(mu[j-1],phi[j-1],x,n)
 	  m.new = m(cand,phi[j-1],x,n)

	  #Draw an observation from a uniform(0,1)

 	  u = runif(1)

	  #If difference in log posterior values is greater than log uniform
          #observation, accept new mu value and increment acceptance counter
	  #by 1.

          if((m.new - m.old) > log(u)) {
		mu[j] = cand
		acceptance.mu = acceptance.mu+1
	  }
  }

  #Metropolis-Hastings step for phi

  #Draw a candidate observation from a normal distribution
  #centered at the previous observation in the chain
  #with a pre-specified variance.

  cand = rnorm(1,phi[j-1],sigma.phi)

  #Check if candidate is between 0 and 1. If not, discard it.
  
  if( (cand > 0) & (cand < 1)) {

	#m.old = Proportional log posterior function at all previous values
	#m.new = Proportional log posterior function at all previous values
	#           except for candidate phi value.
         
 	m.old = m(mu[j-1],phi[j-1],x,n)
	m.new = m(mu[j-1],cand,x,n)


	#Draw an observation from a uniform(0,1)

 	u = runif(1)

	#If difference in log posterior values is greater than log uniform
        #observation, accept new phi value and increment acceptance counter
	#by 1.

        if((m.new - m.old) > log(u)) {
		phi[j] = cand
		acceptance.phi = acceptance.phi + 1
	} 


  }

  #Gibbs sampling step for theta
  #For each player, simulate a batting average from a beta distribution
  #With parameters alpha = x + mu*(1-phi)/phi and beta = n - x + (1-mu)*(1-phi)/phi
  #where mu and phi are the current values in the chain.

  for(i in 1:length(n)) {
   theta[i, j] = rbeta(1, (1-phi[j])/phi[j]*mu[j] + x[i], (1-phi[j])/phi[j]*(1-mu[j]) + n[i] - x[i])
  } 

 } #End MCMC Loop

 #Throw out the burn-in iterations

 mu <- mu[(burn.in + 1):(burn.in + n.draws)]
 phi <- phi[(burn.in + 1):(burn.in + n.draws)]
 theta <- theta[,(burn.in + 1):(burn.in + n.draws)]

 #Return chain for each parameter and acceptance rates for MH step parameters

 return(list(mu = mu, phi = phi, theta = theta, acceptance = c(acceptance.mu/(burn.in + n.draws), acceptance.phi/(burn.in + n.draws))))

}
