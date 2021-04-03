#Data

mu = 0.3306363
M = 302.5670532

x = 114
n = 300

n.new = 250

#Normal-Normal Mean Interval

((n/(n+M))*(x/n) + (M/(n+M))*mu) - 1.96*sqrt(mu*(1-mu)/(M + n))
((n/(n+M))*(x/n) + (M/(n+M))*mu) + 1.96*sqrt(mu*(1-mu)/(M + n))


#Normal Approximation to Beta Mean Interval

(x + mu*M)/(n + M) - 1.96*sqrt((x+mu*M)*(n - x + (1-mu)*M)/((n + M)^2*(1 + n + M)))
(x + mu*M)/(n + M) + 1.96*sqrt((x+mu*M)*(n - x + (1-mu)*M)/((n + M)^2*(1 + n + M)))


#Normal-Normal Prediction Interval

((n/(n+M))*(x/n) + (M/(n+M))*mu) - 1.96*sqrt(mu*(1-mu)/(M + n) + mu*(1-mu)/n.new)
((n/(n+M))*(x/n) + (M/(n+M))*mu) + 1.96*sqrt(mu*(1-mu)/(M + n) + mu*(1-mu)/n.new)

#Normal Approximation to Beta-Binomial Prediction Interval

(x + mu*M)/(n + M) - 1.96*sqrt((x+mu*M)*(n - x + (1-mu)*M)*(n + M + n.new)/(n.new*(n + M)^2*(1 + n + M)))
(x + mu*M)/(n + M) + 1.96*sqrt((x+mu*M)*(n - x + (1-mu)*M)*(n + M + n.new)/(n.new*(n + M)^2*(1 + n + M)))


