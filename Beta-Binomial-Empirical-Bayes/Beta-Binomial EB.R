x <- c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)

y <- round(x*45)
nY <- rep(45, length(y))

theta <- c(.346,.298,.276,.222,.273,.270,.263,.210,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)

#Beta-binomial EB


#Note: when all sample sizes are equal (as here, where nY = 45 fo all observations), then the naive estimator
#below is equivalent to the fancier iterative method I also give

mu = sum(y)/sum(nY)
N = length(y)
s2 = N*sum(nY*(y/nY-mu)^2)/((N-1)*sum(nY))
M = (mu*(1-mu)-s2)/(s2-mu*(1-mu)/N*sum(1/nY))

a = mu * M
b = M - a

B = M/(M + nY)

thetaEB = rep(0, length(y))
for(i in 1:length(y)) thetaEB[i] = (1-B[i])*(y[i]/nY[i]) + B[i]*mu

alphaEB <- y + mu*M
betaEB <- nY - y + (1-mu)*M

#James-Stein Estimator


f = sqrt(45)*asin(2*x-1)
ftheta = sqrt(45)*asin(2*theta - 1)

back = 0.5*sin(f/sqrt(45))+0.5

js = mean(f) + (1 - (18-3)/(sum((f - mean(f))^2)))*(f - mean(f))

backjs = 0.5*sin(js/sqrt(45))+0.5




#Comparison

sum((x - theta)^2)
sum((thetaEB - theta)^2)
sum((backjs - theta)^2)

data.frame("Hits" = y, "MLE" = x, "JS" = backjs, "EB" = thetaEB, "Theta" = theta)








#The code below contains me trying to make confidence intervals based on a procedure
#in a Carl Morris article. It has no relevance to the empirical Bayesian procedure above,
#especially since it uses a normal-normal model, but I thougt I'd leave it in anyway
#just because.

#Intervals

lEB <- qbeta(.025,alphaEB,betaEB)
uEB <- qbeta(.975,alphaEB,betaEB)

fB = (1 - (18-3)/(sum((f - mean(f))^2)))
#fS = sqrt(fB + 3/(18-3)*(f*(1-fB))^2)
fS = sqrt(fB)
ljs = 0.5*sin((js-1.96*fS)/sqrt(45))+0.5
ujs = 0.5*sin((js+1.96*fS)/sqrt(45))+0.5

data.frame(lEB, thetaEB, uEB, ljs, backjs, ujs)


#Morris on transformed

x <- c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)

#x <- .4841 + .0659*sqrt(45)*asin(2*x-1)

f = sqrt(45)*asin(2*x-1)
ftheta = sqrt(45)*asin(2*theta - 1)

back = 0.5*sin(f/sqrt(45))+0.5

js = mean(f) + (1 - (18-3)/(sum((f - mean(f))^2)))*(f - mean(f))

backjs = 0.5*sin(js/sqrt(45))+0.5





S = sum((f - mean(f))^2)
V = 1
m = (18-3)/2
k = 18

BJS = 2*m*V/S
T = S/(2*V)

M = integrate( function(x) exp((1-x)*T)*m*x^(m-1),0,1)$value
B = 2*m*V/S*(1-1/M)

v = 1/m*B^2 - (BJS - B)*(1-(m+1)/m*B)
s2 = V*(1-(k-1)/k*B) + v*(f - mean(f))^2
s = sqrt(s2)


fS = sqrt(fB)
ljs = 0.5*sin((js-1.96*s)/sqrt(45))+0.5
ujs = 0.5*sin((js+1.96*s)/sqrt(45))+0.5

data.frame(lEB, thetaEB, uEB, ljs, backjs, ujs)

s_ = s*sd(x)
ljs = backjs-1.96*s_
ujs = backjs+1.96*s_

data.frame(lEB, thetaEB, uEB, ljs, backjs, ujs)


#Morris on raw data

S = sum((x - mean(x))^2)
V = var(x)
m = (18-3)/2
k = 18

BJS = 2*m*V/S
T = S/(2*V)

M = integrate( function(x) exp((1-x)*T)*m*x^(m-1),0,1)$value
B = 2*m*V/S*(1-1/M)

v = 1/m*B^2 - (BJS - B)*(1-(m+1)/m*B)
s2 = V*(1-(k-1)/k*B) + v*(f - mean(f))^2
s = sqrt(s2)


fS = sqrt(fB)
ljs = x-1.96*s
ujs = x+1.96*s

data.frame(lEB, thetaEB, uEB, ljs, backjs, ujs)

