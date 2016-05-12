#Read in historical data for pitchers

d <- read.csv('C:/Users/Robert/Documents/GitHub/Blog-Code/An-Attempt-To-Create-The-Most-Useless-Baseball-Statistic/Historical BABIP Data.csv',header=T)

#Set start and ending years for calculating year-to-year correlation
#and the minimum number of BIP in a season to be considered in the sample.

startYear = 1901
endYear = 2015
cutoff = 300

#Get the season for each data point

season = d[,1]

#Create empty vectors for the mean and variance of my stat each year,
#the mean and variance of standard normal CDF of normalized BABIP for each
#year (p), and the mean and variance for pitcher BABIP each year (xn)

mu.stat = rep(0, length(startYear:endYear))
var.stat = rep(0, length(startYear:endYear))

mu.p = rep(0, length(startYear:endYear))
var.p = rep(0, length(startYear:endYear))

mu.xn = rep(0, length(startYear:endYear))
var.xn = rep(0, length(startYear:endYear))

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

 p <- 2*pnorm(-abs(((x/n) - mean(x/n))/sd(x/n)))

 #Calculate the average and variance of my new stat 1/(1-p)
 #also calculate the mean and variance of BABIP for that year
 #and the mean and variance of standard normal CDF of normalized
 # BABIP for that year

 mu.stat[year-startYear+1] = mean(1/(1-p))
 var.stat[year-startYear+1] = var(1/(1-p))

 mu.xn[year-startYear+1] = mean(x/n)
 var.xn[year-startYear+1] = var(x/n)

 mu.p[year-startYear+1] = mean(p)
 var.p[year-startYear+1] = var(p)

}

#Create data frame of my statistic
#Tell R to display lots of digits instead of scientific notation

options('scipen' = 999)
data.frame('Year' = startYear:endYear, mu.stat, var.stat)

#Plot mean and variance of my statistic by year

par(mfrow = c(1,2))

plot(startYear:endYear,  mu.stat, type = 'l', xlab = 'Year', ylab = 'Mean', main = 'Mean of Statistic by Year')
plot(startYear:endYear, var.stat, type = 'l', xlab = 'Year', ylab = 'Variance', main = 'Variance of Statistic by Year')



#Year-to-Year correlation


#Define start and end years as before

startYear = 1901
endYear = 2014

#Empty vector of year-to-year correlations

cor.stat <- rep(0, length(startYear:endYear))

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

 babip.1 <- data.frame('Name' = d$Name[samp1], x/n, 'Stat' = 1/(1-p))

 #Calculate statistic for the current year
 #Store the player, BABIP for that year, and stat
 #in a second data frame.

 x <- (d$H - d$HR)[samp2]
 n <- round(((d$H-d$HR)/d$BABIP)[samp2])

 p <- pnorm(((x/n) - mean(x/n))/sd(x/n))

 babip.2 <- data.frame('Name' = d$Name[samp2], x/n, 'Stat' = 1/(1-p))

 #Merge the two into a common data frame by player name

 babip.data = merge(babip.1, babip.2, by = 'Name')

 #Calculate the correlation between this year and the previous year.

 cor.stat[currYear - startYear + 1] = cor(babip.data$Stat.x,babip.data$Stat.y)

}

#Basic R graphics plot of year-to-year correlation of my statistic versus time

data.frame('Year' = startYear:endYear, 'Correlation' = cor.stat)


plot(startYear:endYear,cor.stat, xlab = 'Year')



