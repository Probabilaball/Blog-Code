#Bayesian Plots

#Priors

curve(dbeta(x, 1, 1), xlab = "Batting Average", ylab = "Probability", main = "Prior Belief for Observer A")

curve(dbeta(x, .265*200, (1-.265)*200), xlab = "Batting Average", ylab = "Probability", main = "Prior Belief for Observer B")

#Posteriors

curve(dbeta(x, 1+15, 1+50-15), xlab = "Batting Average", ylab = "Probability", main = "Posterior Belief for Observer A")

curve(dbeta(x, 0.265*200+15, (1-.265)*200+50-15), xlab = "Batting Average", ylab = "Probability", main = "Posterior Belief for Observer B")