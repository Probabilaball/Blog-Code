#Plot 

plot(woba.1, woba.2, pch = 16, cex = 1, xlab = 'First half wOBA', ylab = 'Second half wOBA', ylim = c(0.1, 0.5), col = 'Red')
arrows(woba.1, pred.1$L, woba.1, pred.1$U, length = 0, code=3)