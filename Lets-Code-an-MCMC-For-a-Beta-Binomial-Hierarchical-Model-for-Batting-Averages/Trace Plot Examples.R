#Set the seed for consistent output

set.seed(5)


#Poor mixing


chain.1 <- betaBin.mcmc(x,n, 0.265, 0.00157)
chain.2 <- betaBin.mcmc(x,n, 0.285, 0.00157, sigma.mu = 0.1)
chain.3 <- betaBin.mcmc(x,n, 0.245, 0.00157, sigma.mu = 0.1)


chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")

matplot(data.frame(chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")


#Slow mixing

chain.1 <- betaBin.mcmc(x,n, 0.265, 0.00157)
chain.2 <- betaBin.mcmc(x,n, 0.285, 0.00157, sigma.mu = 0.0001)
chain.3 <- betaBin.mcmc(x,n, 0.245, 0.00157, sigma.mu = 0.0001)

chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")

matplot(data.frame(chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")


#Good mixing

chain.1 <- betaBin.mcmc(x,n, 0.265, 0.00157)
chain.2 <- betaBin.mcmc(x,n, 0.285, 0.00157)
chain.3 <- betaBin.mcmc(x,n, 0.245, 0.00157)

chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")

matplot(data.frame(chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")


