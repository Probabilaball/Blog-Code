#Prior Plots

curve(dbeta(x, 1, 1), xlab = "Batting Average", ylab = "Probability", main = "Prior Belief for Observer A")

curve(dbeta(x, .265*200, (1-.265)*200), xlab = "Batting Average", ylab = "Probability", main = "Prior Belief for Observer B")

#Posterior Plots

curve(dbeta(x, 1+15, 1+50-15), xlab = "Batting Average", ylab = "Probability", main = "Posterior Belief for Observer A")

curve(dbeta(x, 0.265*200+15, (1-.265)*200+50-15), xlab = "Batting Average", ylab = "Probability", main = "Posterior Belief for Observer B")

#Parameters

a = .265*200 + 15
b = (1-.265)*200 + 50-15

#Plot including credible interval

l = qbeta(.025,16, 36)
u = qbeta(.975,16,36)
curve(dbeta(x, 1+15, 1+50-15), xlab = "Batting Average", ylab = "Probability", main = "Posterior Belief for Observer A")
abline(v = l)
abline(v = u)