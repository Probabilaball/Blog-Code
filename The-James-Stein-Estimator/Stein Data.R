x = c(.400,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)

theta = c(.346,.298,.276,.222,.273,.270,.263,.210,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)

f = sqrt(45)*asin(2*x-1)
ftheta = sqrt(45)*asin(2*theta - 1)

back = 0.5*sin(f/sqrt(45))+0.5

js = mean(f) + (1 - (18-3)/(sum((f - mean(f))^2)))*(f - mean(f))

sum((f - ftheta)^2)
sum((js - ftheta)^2)


backjs = 0.5*sin(js/sqrt(45))+0.5

sum((x - theta)^2)
sum((backjs - theta)^2)

data.frame(round(x,3), round(f,3), round(js,3), round(backjs,3), theta)