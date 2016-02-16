
theta <- c(.346,.298,.276,.222,.273,.270,.263,.210,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)
x <- c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)
x <- round(x*45)
n <- rep(45, 18)


par <- mm.est(x,n)

res <- mm.bootstrap(x,n)


mean((res$Naive[,1] < theta) & (theta < res$Naive[,2]))

mean((res$Bootstrap[,1] < theta) & (theta < res$Bootstrap[,2]))


#Plots of player 12


 alpha = par$alpha
 beta = par$beta

 alphaBoot <- rep(0, 5000)
 betaBoot <- rep(0, 5000)

 for(i in 1:5000) {

  #Generate data from two-stage model
  #Using parameter estimates from method of moments estimator

  theta <- rbeta(length(x), alpha, beta)
  xBoot <- rbinom(length(x), n, theta)

  #Store bootstrapped parameter estimates in appropriate vectors

  parBoot <- mm.est(xBoot, n)

  alphaBoot[i] <- parBoot$alpha
  betaBoot[i] <- parBoot$beta

 }

 xPts <- seq(0.001,0.999,by=0.001)
 yPts <- rep(0, length(xPts))

 for(i in 1:length(xPts)) yPts[i] = mean(dbeta(xPts[i], 10 + alphaBoot, 45-10 + betaBoot),na.rm=T)

 
 curve(dbeta(x, 10 + par$alpha, 45 - 10 + par$beta), xlim = c(0.15, 0.4), xlab = 'True Batting Average', ylab = 'Posterior Density', main = "Naive and Bootstrap Corrected Posterior Densities\nFor 10 Hits in 45 At-Bats")
 points(xPts, yPts, type = 'l',lty = 2)
