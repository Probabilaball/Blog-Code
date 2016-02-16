#Load Data

d <- read.csv('C:/[File Path]/Interval Data.csv',header=T)

#Interval Calculation Function
#x.1 and n.1 are the count of events and trials in the first half
#x.2 and n.2 are the count of events and trials in the second half
#mu is the population mean
#M is the stabilization point

calcInt <- function(x.1, n.1, x.2, n.2, mu, M) {

 #Get the total number of observations

 N <- length(x.1)


 #Normal-Normal Mean Intervals


 normal.mean.L <- ((n.1/(M+n.1))*(x.1/n.1) + (M/(M+n.1))*mu)-1.96*sqrt(mu*(1-mu)/(M+n.1))
 normal.mean.U <- ((n.1/(M+n.1))*(x.1/n.1) + (M/(M+n.1))*mu)+1.96*sqrt(mu*(1-mu)/(M+n.1))


 normal.mean.coverage <- mean((normal.mean.L <= x.2/n.2) & (x.2/n.2 <= normal.mean.U))
 normal.mean.length <- mean(normal.mean.U-normal.mean.L)


 #Normal-Normal Prediction Intervals

 normal.pred.L <- ((n.1/(M+n.1))*(x.1/n.1) + (M/(M+n.1))*mu)-1.96*sqrt(mu*(1-mu)/(M+n.1) + mu*(1-mu)/n.2)
 normal.pred.U <- ((n.1/(M+n.1))*(x.1/n.1) + (M/(M+n.1))*mu)+1.96*sqrt(mu*(1-mu)/(M+n.1) + mu*(1-mu)/n.2)

 normal.pred.coverage <- mean((normal.pred.L <= x.2/n.2) & (x.2/n.2 <= normal.pred.U))
 normal.pred.length <- mean(normal.pred.U-normal.pred.L)


 #Beta-Binomial Mean Intervals


 betabin.mean.L <- rep(0, N)
 betabin.mean.U <- rep(0, N)

 for(i in 1:N) {

  int <- qbeta(c(.025,.975), x.1[i] + mu*M, n.1[i]-x.1[i] + (1-mu)*M)

  betabin.mean.L[i] <- int[1]
  betabin.mean.U[i] <- int[2]

 }

 betabin.mean.coverage <- mean((betabin.mean.L <= x.2/n.2) & (x.2/n.2 <= betabin.mean.U))
 betabin.mean.length <- mean(betabin.mean.U-betabin.mean.L)


 #Beta-Binomial Mean Intervals (Normal Approximation)

 betabin.norm.mean.L <- (x.1 + mu*M)/(n.1 + M)-1.96*sqrt(((x.1+mu*M)*(n.1-x.1+(1-mu)*M))/((n.1+M)^2*(1+n.1+M)))
 betabin.norm.mean.U <- (x.1 + mu*M)/(n.1 + M)+1.96*sqrt(((x.1+mu*M)*(n.1-x.1+(1-mu)*M))/((n.1+M)^2*(1+n.1+M)))

 betabin.norm.mean.coverage <- mean((betabin.norm.mean.L <= x.2/n.2) & (x.2/n.2 <= betabin.norm.mean.U))
 betabin.norm.mean.length <- mean(betabin.norm.mean.U-betabin.norm.mean.L)


 #Beta-Binomial Prediction Intervals

 #Beta-binomial density function

 dbetabinom <- function(x, n, alpha, beta) choose(n,x)*beta(x+alpha, n-x+beta)/beta(alpha, beta)

 #Beta-binomial quantile function
 #Simple search - start at 0 and keep count of cumulative probability
 #until you reach the desired quantile.
 #There is probably a better way to do this.

 qbetabinom <- function(q,n, alpha, beta) {
  
  x <- -1
  p <- 0
 
  while(p < q) {
 
   x <- x + 1
   p <- p + dbetabinom(x, n, alpha, beta)

  }

  return(x)

 }


 betabin.pred.L <- rep(0, N)
 betabin.pred.U <- rep(0, N)

 for(i in 1:N) {

  int <- c(0,0)
  int[1] <- qbetabinom(.025, n.2[i], x.1[i] + mu*M, n.1[i]-x.1[i] + (1-mu)*M)
  int[2] <- qbetabinom(.975, n.2[i], x.1[i] + mu*M, n.1[i]-x.1[i] + (1-mu)*M)

  betabin.pred.L[i] <- int[1]/n.2[i]
  betabin.pred.U[i] <- int[2]/n.2[i]

 }

 betabin.pred.coverage <- mean((betabin.pred.L <= x.2/n.2) & (x.2/n.2 <= betabin.pred.U))
 betabin.pred.length <- mean(betabin.pred.U-betabin.pred.L)


 #Beta-Binomial Prediction Intervals (Normal Approximation)

 betabin.norm.pred.L <- (x.1 + mu*M)/(n.1 + M)-1.96*sqrt(((x.1+mu*M)*(n.1-x.1+(1-mu)*M)*(n.1 + M + n.2))/(n.2*(n.1+M)^2*(1+n.1+M)))
 betabin.norm.pred.U <- (x.1 + mu*M)/(n.1 + M)+1.96*sqrt(((x.1+mu*M)*(n.1-x.1+(1-mu)*M)*(n.1 + M + n.2))/(n.2*(n.1+M)^2*(1+n.1+M)))

 betabin.norm.pred.coverage <- mean((betabin.norm.pred.L <= x.2/n.2) & (x.2/n.2 <= betabin.norm.pred.U))
 betabin.norm.pred.length <- mean(betabin.norm.pred.U-betabin.norm.pred.L)

 #Return Results
 #List of lists
 #Each list contains four elements: lower interval bounds, upper interval bounds
 #coverage, and average length

 normal.mean <- list('lower' = normal.mean.L, 'upper' = normal.mean.U, 'coverage' = normal.mean.coverage, 'length' = normal.mean.length)
 betabin.mean <- list('lower' = betabin.mean.L, 'upper' = betabin.mean.U, 'coverage' = betabin.mean.coverage, 'length' = betabin.mean.length)
 betabin.norm.mean <- list('lower' = betabin.norm.mean.L, 'upper' = betabin.norm.mean.U, 'coverage' = betabin.norm.mean.coverage, 'length' = betabin.norm.mean.length)

 normal.pred <- list('lower' = normal.pred.L, 'upper' = normal.pred.U, 'coverage' = normal.pred.coverage, 'length' = normal.pred.length)
 betabin.pred <- list('lower' = betabin.pred.L, 'upper' = betabin.pred.U, 'coverage' = normal.pred.coverage, 'length' = normal.pred.length)
 betabin.norm.pred <- list('lower' = betabin.norm.pred.L, 'upper' = betabin.norm.pred.U, 'coverage' = betabin.norm.pred.coverage, 'length' = betabin.norm.pred.length)


 return(list('normal.mean' = normal.mean, 'betabin.mean' = betabin.mean, 'betabin.norm.mean' = betabin.norm.mean, 'normal.pred' = normal.pred, 'betabin.pred' = betabin.pred, 'betabin.norm.pred' = betabin.norm.pred))

}

#On-Base Percentage Calculations

x.1 <- d$H.1 + d$BB.1 + d$HBP.1
n.1 <- d$PA.1

x.2 <- d$H.2 + d$BB.2 + d$HBP.2
n.2 <- d$PA.2

muOBP <- 0.330
MOBP <- 296

resOBP <- calcInt(x.1, n.1, x.2, n.2, muOBP, MOBP)

#Batting Average Calculations


x.1 <- d$H.1
n.1 <- d$AB.1

x.2 <- d$H.2
n.2 <- d$AB.2

muBA <- 0.268
MBA <- 466

resBA <- calcInt(x.1, n.1, x.2, n.2, muBA, MBA)

#1B Rate Calculations

x.1 <- d$X1B.1
n.1 <- d$PA.1

x.2 <- d$X1B.2
n.2 <- d$PA.2

mu1B <- 0.158
M1B <- 222

res1B <- calcInt(x.1, n.1, x.2, n.2, mu1B, M1B)

#2B Rate Calculations

x.1 <- d$X2B.1
n.1 <- d$PA.1

x.2 <- d$X2B.2
n.2 <- d$PA.2

mu2B <- .0475
M2B <- 1025

res2B <- calcInt(x.1, n.1, x.2, n.2, mu2B, M2B)

#3B Rate Calculations

x.1 <- d$X3B.1
n.1 <- d$PA.1

x.2 <- d$X3B.2
n.2 <- d$PA.2

mu3B <- 0.00492
M3B <- 373

res3B <- calcInt(x.1, n.1, x.2, n.2, mu3B, M3B)

#XBH Rate Calculations

x.1 <- d$X3B.1 + d$X2B.1
n.1 <- d$PA.1

x.2 <- d$X3B.2 + d$X2B.2
n.2 <- d$PA.2

muXBH <- 0.0524
MXBH <- 1006

resXBH <- calcInt(x.1, n.1, x.2, n.2, muXBH, MXBH)

#HR Rate Calculations

x.1 <- d$HR.1
n.1 <- d$PA.1

x.2 <- d$HR.2
n.2 <- d$PA.2

muHR <- .0274
MHR <- 125

resHR <- calcInt(x.1, n.1, x.2, n.2, muHR, MHR)

#SO Rate Calculations

x.1 <- d$SO.1
n.1 <- d$PA.1

x.2 <- d$SO.2
n.2 <- d$PA.2

muSO <- 0.181
MSO <- 50

resSO <- calcInt(x.1, n.1, x.2, n.2, muSO, MSO)

#BB Rate Calculations

x.1 <- d$BB.1
n.1 <- d$PA.1

x.2 <- d$BB.2
n.2 <- d$PA.2

muBB <- 0.085
MBB <- 106

resBB <- calcInt(x.1, n.1, x.2, n.2, muBB, MBB)

#HBP Rate Calculations

x.1 <- d$HBP.1
n.1 <- d$PA.1

x.2 <- d$HBP.2
n.2 <- d$PA.2

muHBP <- 0.00866
MHBP <- 297

resHBP <- calcInt(x.1, n.1, x.2, n.2, muHBP, MHBP)

#Create Table code for blog

#Mean Intervals

stats <- c('OBP', 'BA', '1B', '2B', '3B', 'XBH', 'HR','BB','SO','HBP')

for(i in 1:length(stats)) {

 var.res <- paste('res', stats[i], sep='')
 
 mu <- get(paste('mu', stats[i], sep = ''))
 M <- get(paste('M', stats[i], sep = ''))

 s <- paste(stats[i], mu, M, round(get(var.res)[[1]][[3]],3),round(get(var.res)[[2]][[3]],3),round(get(var.res)[[3]][[3]],3), sep = ' & ')
 s <- paste(s, '\\\\')
 cat(s, "\n")

} 


#Prediction Intervals

stats <- c('OBP', 'BA', '1B', '2B', '3B', 'XBH', 'HR','BB','SO','HBP')

for(i in 1:length(stats)) {

 var.res <- paste('res', stats[i], sep='')
 
 mu <- get(paste('mu', stats[i], sep = ''))
 M <- get(paste('M', stats[i], sep = ''))

 s <- paste(stats[i], mu, M, round(get(var.res)[[4]][[3]],3),round(get(var.res)[[5]][[3]],3),round(get(var.res)[[6]][[3]],3), sep = ' & ')
 s <- paste(s, '\\\\')
 cat(s, "\n")

} 
