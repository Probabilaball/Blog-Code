#Read in historical data for pitchers

d <- read.csv('C:/Users/Robert/Dropbox/Baseball Blog Articles/An Attempt to Create the Most Useless Baseball Stat/Historical BABIP Data.csv',header=T)

#Set start and ending years for calculating year-to-year correlation
#and the minimum number of BIP in a season to be considered in the sample.

startYear = 1901
endYear = 2014
cutoff = 300

#Get the season for each data point

season = d[,1]

#Create empty vectors for the mean and variance of my stat each year (LIFT),
#the mean and variance of standard normal CDF of normalized BABIP for each
#year (P), and the mean and variance for pitcher BABIP each year (XN)

muLIFT = rep(0, length(startYear:endYear))
varLIFT = rep(0, length(startYear:endYear))

muP = rep(0, length(startYear:endYear))
varP = rep(0, length(startYear:endYear))

muXN = rep(0, length(startYear:endYear))
varXN = rep(0, length(startYear:endYear))

#Begin calculation loop

for(year in startYear:endYear) {

 #Only take the current year

 samp = (season == year)

 #Get number BIP by rounding and H - HR by exact calculation

 x <- (d$H - d$HR)[samp]
 n <- round(((d$H-d$HR)/d$BABIP)[samp])

 #Get rid of all the NAs and pitchers that didn't have enough BIP

 s = (!is.na(x)) & (!is.na(n)) & (n >= cutoff)

 x <- x[s]
 n <- n[s]

 #Standardize the BABIP for each pitcher by subtracting the mean BABIP
 #for that year and dividing by the standard deviation of BABIP for that year
 #then take the standard normal CDF of those values 

 p <- pnorm(((x/n) - mean(x/n))/sd(x/n))

 #Calculate the average and variance of my new stat 1/(1-p)
 #also calculate the mean and variance of BABIP for that year
 #and the mean and variance of standard normal CDF of normalized
 # BABIP for that year

 muLIFT[year-startYear+1] = mean(1/(1-p))
 varLIFT[year-startYear+1] = var(1/(1-p))

 muXN[year-startYear+1] = mean(x/n)
 varXN[year-startYear+1] = var(x/n)

 muP[year-startYear+1] = mean(p)
 varP[year-startYear+1] = var(p)

}

#Create data frame of my statistic
#Tell R to display lots of digits instead of scientific notation


options('scipen' = 999)
data.frame('Year' = startYear:endYear, muLIFT, varLIFT)

#Plot mean and variance of my statistic by year

plot(startYear:endYear, muLIFT, type = 'l')
plot(startYear:endYear, varLIFT, type = 'l')



#Year-to-Year correlation


#Define start and end years as before

startYear = 1901
endYear = 2014

#Empty vector of year-to-year correlations

corLIFT <- rep(0, length(startYear:endYear))

#Minimum number of BIP to be included in sample

cutoff <- 300


for(currYear in startYear:endYear) {

 #Take samples of players in current year and previous year who
 #meet cutoff requirements and don't have any NAs in their data

 samp1 <- (season == (currYear - 1)) & ((d$H - d$HR)/d$BABIP >= cutoff) & !is.na((d$H - d$HR)/d$BABIP)
 samp2 <- (season == currYear) & ((d$H-d$HR)/d$BABIP >= cutoff)& !is.na((d$H - d$HR)/d$BABIP)


 #Calculate statistic for the previous year
 #Store the player, BABIP for that year, and stat
 #in a data frame.

 x <- (d$H - d$HR)[samp1]
 n <- round(((d$H-d$HR)/d$BABIP)[samp1])

 p <- pnorm(((x/n) - mean(x/n))/sd(x/n))

 babip.1 <- data.frame('Name' = d$Name[samp1], x/n, 'LIFT' = 1/(1-p))

 #Calculate statistic for the current year
 #Store the player, BABIP for that year, and stat
 #in a second data frame.

 x <- (d$H - d$HR)[samp2]
 n <- round(((d$H-d$HR)/d$BABIP)[samp2])

 p <- pnorm(((x/n) - mean(x/n))/sd(x/n))

 babip.2 <- data.frame('Name' = d$Name[samp2], x/n, 'LIFT' = 1/(1-p))

 #Merge the two into a common data frame by player name

 babip.data = merge(babip.1, babip.2, by = 'Name')

 #Calculate the correlation between this year and the previous year.

 corLIFT[currYear - startYear + 1] = cor(babip.data$LIFT.x,babip.data$LIFT.y)

}

#Basic R graphics plot of year-to-year correlation of my statistic versus time

plot(startYear:endYear,corLIFT, xlab = 'Year')

#Fancier ggplot of year-to-year correlation versus time

library(ggplot2)
qplot(startYear:endYear, corLIFT,xlab = 'Year', ylab = 'Year-to-year Correlation')

