f = function(x) x*(1-x)^2

x = seq(0,1,by=0.001)
y = rep(0, length(x))

for(i in 1:length(x)) y[i] = f(x[i])

plot(x,y, type = 'l')
