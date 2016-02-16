MSO = 49.73072
MBB = 105.5944
MHR = 124.5226
M1B = 222.1613
MOBP = 295.7912
MHBP = 297.4059
MXBH = 1006.301


#Based on p

p <- seq(0.4,.8,by=0.001)
p_ <- p/(1-p)

sPoints <- cbind(p_*MSO, p_*MBB, p_*MHR, p_*M1B, p_*MOBP, p_*MHBP,p_*MXBH)

matplot(p, sPoints, type = 'l', xlab = "Stabilization Level", ylab = "Sample Size Required", col=c(1:6,8), lty = 1:7)

legend(0.4, 4000,  c("SO Rate", "BB Rate", "HR Rate", "1B Rate", "OBP Rate", "HBP Rate", "XBH Rate"), col=c(1:6,8), lty = 1:7)


#Based on n

n <- 1:600

#sPoints <- cbind(n/(n+MOBP), n/(n+MSO), n/(n+MBB), n/(n+M1B), n/(n+MXBH), n/(n+MHR))
sPoints <- cbind(n/(n+MSO), n/(n+MBB), n/(n+MHR), n/(n+M1B), n/(n+MOBP), n/(n+MHBP),n/(n+MXBH))


matplot(n, sPoints, type = 'l', xlab = "Sample size", ylab = "Stabilization Level", col=c(1:6,8),lty=1:7)

legend(400, 0.37, c("SO Rate", "BB Rate", "HR Rate", "1B Rate", "OBP Rate", "HBP Rate","XBH Rate"), col = c(1:6,8), lty = 1:7)

#Order: SO Rate, BB Rate, HR Rate, 1B Rate, OBP , HBP, XBH