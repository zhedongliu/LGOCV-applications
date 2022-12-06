n = 1e5
T1 = rweibull(n,shape = 1,scale = 1)
T2 = rweibull(n,shape = 2,scale = 1)
T3 = rweibull(n,shape = 2,scale = 1)

Time = numeric(n)
Event = numeric(n)
for(i in 1:n){
    Time[i] = min(c(T1[i],T2[i],T3[i]))
    Event[i] = which.min(c(T1[i],T2[i],T3[i]))
}

tmax = max(T1,T2,T3)
tmin = min(T1,T2,T3)

t = seq(0,tmax,0.01)
