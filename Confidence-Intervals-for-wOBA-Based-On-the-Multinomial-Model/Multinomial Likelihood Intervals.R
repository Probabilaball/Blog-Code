ml.multinomial <- function(data, n) {

 par <- data/n

 k <- length(data)

 v <- matrix(rep(0, k*k),k,k)

 for(i in 1:k) {
  for(j in 1:k) { 
   v[i,j] = -par[i]*par[j]/n
  }
 }

 for(i in 1:k) v[i,i] = par[i]*(1-par[i])/n

 return(list(par = par, v = v))


}

#Mike Trout
#Data = 2013 Season

uBB = 100 #unintentional walks
HBP = 9# HBP
S = 115# singles
D = 39# doubles
T = 9# triples
HR = 27# home runs

n = 716-10 #PA - IBB

res <- ml.multinomial(c(uBB, HBP, S, D, T, HR), n)
par <- res$par
v <- res$v

#wOBA = (0.690*uBB + 0.722*HBP + 0.888*S + 1.271*D + 1.616*T +2.101*HR)/n
w = matrix(c(.69, .722, .888, 1.271, 1.616, 2.101),1,6)

wOBA = sum(w*par)
wOBA


#Confidence Interval

wOBA - 1.96*sqrt(w%*%v%*%t(w))
wOBA + 1.96*sqrt(w%*%v%*%t(w))


