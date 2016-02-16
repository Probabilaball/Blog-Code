h =100
n = 400
p = h/n

l = p - 1.96*sqrt(p*(1-p)/n)
u = p + 1.96*sqrt(p*(1-p)/n)

o = p/(1-p)
vo = p/(n*(1-p)^3)

lo = o - 1.96*sqrt(vo)
uo = o + 1.96*sqrt(vo)