#BIP Stats

MFB = 61.95882
MGB = 65.51318
MLD = 768.4226
MBABIP = 2006.711

#Based on n

n <- 1:1000

sPoints <- cbind(n/(n+MFB), n/(n+MGB), n/(n+MLD), n/(n+MBABIP))


matplot(n, sPoints, type = 'l', xlab = "Balls in Play", ylab = "Stabilization Level", col=c(1:6,8),lty=1:7)

legend(700, 0.85, c("FB Rate", "GB Rate", "LD Rate", "BABIP"), col = c(1:6,8), lty = 1:7)


#TBF Stats

MSO = 90.93711
MBB = 221.2492
MOBP = 524.731
MH = 623.3454
MHR = 931.5889
MHBP = 989.304



#Based on n

n <- 1:1000

sPoints <- cbind(n/(n+MSO), n/(n+MBB), n/(n+MOBP), n/(n+MH), n/(n+MHR), n/(n+MHBP))


matplot(n, sPoints, type = 'l', xlab = "Total Batters Faced", ylab = "Stabilization Level", col=c(1:6,8),lty=1:7)

legend(700, 0.4, c("SO Rate", "BB Rate", "OBP", "Hit Rate", "HR Rate", "HBP Rate" ), col = c(1:6,8), lty = 1:7)
