n_person = 10
age = sample(x = 15:75,size = n_person,replace = T)
age = scale(age)
ui = rnorm(n_person)
eta_event1 = ui
eta_event2 = ui
eta_event3 = ui
scale1 = exp(-eta_event1)
scale2 = exp(-eta_event2)
scale3 = exp(-eta_event3)
T1 = rweibull(n = n_person,shape = 3,scale = scale1)
T2 = rweibull(n = n_person,shape = 3,scale = scale2)
T3 = rweibull(n = n_person,shape = 3,scale = scale3)
Time = numeric(n_person)
Event = numeric(n_person)
for(i in 1:n_person){
    Time[i] = min(T1[i],T2[i],T3[i])
    Event[i] = which.min(c(T1[i],T2[i],T3[i]))
}

s = 0
t = 0.5 #prediction time
pi1 = function(t,scale) dweibull(x = t,shape = 3,scale = scale[1])*pweibull(q = t,shape = 3,scale = scale[2],lower.tail = FALSE)*pweibull(q = t,shape = 3,scale = scale[3],lower.tail = FALSE)
pi2 = function(t,scale) dweibull(x = t,shape = 3,scale = scale[2])*pweibull(q = t,shape = 3,scale = scale[1],lower.tail = FALSE)*pweibull(q = t,shape = 3,scale = scale[3],lower.tail = FALSE)
pi3 = function(t,scale) dweibull(x = t,shape = 3,scale = scale[3])*pweibull(q = t,shape = 3,scale = scale[1],lower.tail = FALSE)*pweibull(q = t,shape = 3,scale = scale[2],lower.tail = FALSE)


p0 = function(t,scale) pweibull(q = t,shape = 3,scale = scale[1],lower.tail = FALSE)*pweibull(q = t,shape = 3,scale = scale[2],lower.tail = FALSE)*pweibull(q = t,shape = 3,scale = scale[3],lower.tail = FALSE)
p1 = function(s,t,scale) quad(f = pi1,xa = s,xb = t,scale = scale)
p2 = function(s,t,scale) quad(f = pi2,xa = s,xb = t,scale = scale)
p3 = function(s,t,scale) quad(f = pi3,xa = s,xb = t,scale = scale)

p = matrix(0,n_person,4)
for(i in 1:n_person){
    p[i,1] = p1(s,t,c(scale1[i],scale2[i],scale3[i]))
    p[i,2] = p2(s,t,c(scale1[i],scale2[i],scale3[i]))
    p[i,3] = p3(s,t,c(scale1[i],scale2[i],scale3[i]))
    p[i,4] = p0(t,c(scale1[i],scale2[i],scale3[i]))/p0(s,c(scale1[i],scale2[i],scale3[i]))
    #print(i)
}


AUC = numeric(4)
for(i in 1:3){
    true_fold = which(Time <= t & Time >= s & Event == i)
    false_fold = which(Time > t | (Time >= s & Event != i))
    for(j in true_fold){
        AUC[i] = AUC[i] + sum(p[j,i] > p[false_fold,i])
    }
    AUC[i] = AUC[i]/length(true_fold)/length(false_fold)
}
true_fold = which(Time > t)
false_fold = which(Time <= t & Time > s)
for(j in true_fold){
    AUC[4] = AUC[4] + sum(p[j,4] > p[false_fold,4])
}
AUC[4] = AUC[4]/length(true_fold)/length(false_fold)
AUC





TPR = function(c,i){
    sum((p[,i]>c)*(Time <= t & Time >= s & Event == i))/sum(Time <= t & Time >= s & Event == i)
}
FPR = function(c,i){
    sum((p[,i]>c)*(Time > t | (Time >= s & Event != i)))/sum(Time > t | (Time >= s & Event != i))
}
c = seq(1,0,-0.00001)
xx = unlist(lapply(c,FPR,i = 1))
yy = unlist(lapply(c,TPR,i=1))
plot(xx,yy,type="l")
pracma::trapz(xx,yy)








