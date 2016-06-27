#Read in Data and separate into two data frames
#Representing the first and second halves of the season

d <- read.csv('C:/Users/Robert/Desktop/Multinomial-Dirichlet Empirical Bayes/Interval Data.csv',header=T)

x.1 <- data.frame(d$X1B.1, d$X2B.1, d$X3B.1, d$HR.1)
x.1 <- as.matrix(cbind(x.1, d$AB.1 - apply(x.1, 1, sum)))

x.2 <- data.frame(d$X1B.2, d$X2B.2, d$X3B.2, d$HR.2)
x.2 <- as.matrix(cbind(x.2, d$AB.2 - apply(x.2, 1, sum)))


#Define variance approximation functions

#Covariance matrix for normal approximation
#to posterior interval for wOBA

cov.matrix <- function(pars) {
  
   k <- length(pars)
   m <- matrix(rep(0,k*k),k,k)

   a0 <- sum(pars)

   for(i in 1:k) {
    for(j in 1:k) {
      if(i == j) m[i,j] = pars[i]*(a0-pars[i])/(a0^2*(a0 + 1))
      if(i != j) m[i,j] = -pars[i]*pars[j]/(a0^2*(a0 + 1))
    }
   }

   return(m)
}

#Covariance matrix for normal approximation
#to posterior predictive interval for wOBA
#in n.new events (PA - IBB)

cov.pred.matrix <- function(pars,n.new) {
  
   k <- length(pars)
   m <- matrix(rep(0,k*k),k,k)

   a0 <- sum(pars)

   for(i in 1:k) {
    for(j in 1:k) {
      if(i == j) m[i,j] = 1/n.new*pars[i]/a0*(1-pars[i]/a0)*(n.new+a0)/(1 + a0)
      if(i != j) m[i,j] = -1/n.new*pars[i]*pars[j]/a0^2*(n.new+a0)/(1+a0)
    }
   }

   return(m)
}


#Calculate approximate posterior intervals for the dirichlet-multinomial
#model with data x, weights w, and dirichlet prior parameters par

dirmult.mean.int <- function(x, w, par) {

  N <- dim(x)[1] 

  n <- apply(x, 1, sum)

  stat <- rep(0, N)
  stat.L <- rep(0, N)
  stat.U <- rep(0, N)

  for(i in 1:N) {

   t <- x[i,] + par

   stat[i] <- sum(w*t)/sum(t)
   
   v <- t(w)%*%cov.matrix(t)%*%w
   
   stat.L[i] <- stat[i] - 1.96*sqrt(v)
   stat.U[i] <- stat[i] + 1.96*sqrt(v)

  }

  return(list('Mean' = stat, 'Lower' = stat.L, 'Upper' = stat.U))

}
   
#Calculate approximate posterior predictive intervals for the 
#dirichlet-multinomial model with data x, weights w, and dirichlet
#prior parameters par, in n.new events

dirmult.pred.int <- function(x, w, par, n.new) {

  N <- dim(x)[1] 

  n <- apply(x, 1, sum)

  stat <- rep(0, N)
  stat.L <- rep(0, N)
  stat.U <- rep(0, N)

  for(i in 1:N) {

   t <- x[i,] + par

   stat[i] <- sum(w*t)/sum(t)
   
   v <- t(w)%*%cov.pred.matrix(t, n.new[i])%*%w
   
   stat.L[i] <- stat[i] - 1.96*sqrt(v)
   stat.U[i] <- stat[i] + 1.96*sqrt(v)

  }

  return(list('Mean' = stat, 'Lower' = stat.L, 'Upper' = stat.U))

}
   
#Calculate approximate intervals for each player and estimate coverage

#Use parameters estimated in previous article from
#2010 - 2015 data with at least 300 events

par <- c(42.443604, 12.855782, 1.381905, 7.073672, 176.120837)

#wOBA weights

w <- c(1,2,3,4,0)

#n.2 = second half events. I'm using this as n.new in the predictive intervals

n.2 <- d$AB.2

#Calculate mean and predictive wOBA intervals for the second half
#Using only the hitting stats in the first half

int.1 <- dirmult.mean.int(x.1, w, par)
pred.1 <- dirmult.pred.int(x.1, w, par, n.2)

#Calculate second-half SLG

slg.2 <- rep(0, dim(x.1)[1])
for(i in 1:length(slg.2)) slg.2[i] <- sum(w*x.2[i,]/n.2[i])

#Calculate coverage of the mean and predictive SLG intervals

mean( (int.1$L <= slg.2) & (slg.2 <= int.1$U) )
mean( (pred.1$L <= slg.2) & (slg.2 <= pred.1$U) )

  