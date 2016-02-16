#Analytic Posterior Predictive

f = function(y) choose(250,y)*beta(76 + y, 426-y)/beta(76,176)

y <- 40:110
p <- rep(0, length(y))

for(i in 1:length(y)) p[i] = f(y[i])

plot(y/250,p, type = 'h', xlab = 'New OBP in 250 PA', ylab = "Probability", main = "Posterior Predictive Distribution for OBP in 250 PA")

#Simulated Posterior Predictive

theta <- rbeta(1000000, 76,176)
y <- rbinom(1000000, 250, theta)

plot(table(y)/(250*length(y)),xlim=c(40,110), xlab = 'New OBP in 250 PA', ylab = "Simulated Probability", main = "Simulated Posterior Predictive Distribution for OBP in 250 PA") 

