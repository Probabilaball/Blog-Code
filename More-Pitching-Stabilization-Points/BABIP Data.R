data <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/More Pitcher Stabilization Points/Pitcher Data.csv',header=T)

x <- data$H - data$HR
n <- data$GB + data$FB + data$LD
g <- data$BABIP


samp <- (n >= 300)

x <- x[samp]
n <- n[samp]
g <- g[samp]

for(i in 1:length(x)) {

  while(abs(round((x/n)[i],3) != g[i]) ) {
   if(round((x/n)[i],3) < g[i]) n[i] = n[i]-1
   if(round((x/n)[i],3) > g[i]) n[i] = n[i]+1
  }
}


muStart = sum(x)/sum(n)
N = length(x)
s2 = N*sum(n*(x/n-muStart)^2)/((N-1)*sum(n))
MStart = (muStart*(1-muStart)-s2)/(s2-muStart*(1-muStart)/N*sum(1/n))
phiStart = 1/(MStart+1)

#Maximize it

ml = mlBetaBinom(c(muStart,phiStart),x,n)
mu = ml$par[1]
phi = ml$par[2]
M = (1-phi)/phi

yMax <- dbeta((mu*M-1)/(M-2), mu*M, (1-mu)*M)

hist(x/n,freq=F, ylim = c(0, yMax), xlab = "BABIP",ylab = "Density", main = "Observed BABIP with True Talent Distribution")
curve(dbeta(x, mu*M, (1-mu)*M), add=T, lty = 2)

v <- (-1/phi^2)^2*ml$v[2,2]

M

sqrt(v)