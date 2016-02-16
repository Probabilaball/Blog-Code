x = c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)

theta = c(.346,.298,.276,.222,.273,.270,.263,.210,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)

#Normal-Normal EB

B = (mean(x)*(1-mean(x))/45)/var(x)

thetaNN = mean(x) + (1-B)*(x-mean(x))

#Beta-binomial EB

y <- round(x*45)
nY <- rep(45, length(y))

mu = sum(y)/sum(nY)
N = length(y)
s2 = N*sum(nY*(y/nY-mu)^2)/((N-1)*sum(nY))
M = (mu*(1-mu)-s2)/(s2-mu*(1-mu)/N*sum(1/nY))

a = mu * M
b = M - a

B = M/(M + nY)

thetaBB = rep(0, length(y))
for(i in 1:length(y)) thetaBB[i] = (1-B[i])*(y[i]/nY[i]) + B[i]*mu

alphaEB <- y + mu*M
betaEB <- nY - y + (1-mu)*M

#James-Stein Estimator


f = sqrt(45)*asin(2*x-1)
ftheta = sqrt(45)*asin(2*theta - 1)

back = 0.5*sin(f/sqrt(45))+0.5

js = mean(f) + (1 - (18-3)/(sum((f - mean(f))^2)))*(f - mean(f))

backjs = 0.5*sin(js/sqrt(45))+0.5

#Comparisons

sum((x - theta)^2)
sum((thetaNN - theta)^2)
sum((thetaBB - theta)^2)
sum((backjs - theta)^2)


data.frame(round(x,3), round(thetaNN,3), round(thetaBB,3), round(backjs,3),theta)