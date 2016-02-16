n = 50

#data <- rbinom(1, n, .3)
data <- 15

l <- function(p, n, data) data*log(p) + (n-data)*log(1-p)

x <- seq(.0001,.9999, by = 0.0001) 
y <- rep(0,length(x))

for(i in 1:length(x)) y[i] = l(x[i], n, data)

plot(x,y, type = 'l', xlab = "Batting Average", ylab = "Log Likelihood")


plot(x,-2*(y - max(y)), xlim = c(data/n-.15, data/n+.15), type = 'l', ylim = c(0,4), xlab = "Theta0", ylab = "Delta Function")
abline(h = qchisq(.95,1), col = 'blue')

min(x[2*(max(y) - y) <= 3.84])
max(x[2*(max(y) - y) <= 3.84])

p = data/n
p - 1.96*sqrt(p*(1-p)/n)
p + 1.96*sqrt(p*(1-p)/n)