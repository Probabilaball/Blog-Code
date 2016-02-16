curve(dbeta(x, 6+36.5,10-6+36.5), xlab = "True Winning Proportion", ylab = "Probability", main = "Posterior of Winning Proportion for a 6-4 Team")

abline(v = (6+36.5)/(10+73),lty = 2)
abline(v = qbeta(.025,6+36.5,10-6+36.5))
abline(v = qbeta(.975,6+36.5,10-6+36.5))