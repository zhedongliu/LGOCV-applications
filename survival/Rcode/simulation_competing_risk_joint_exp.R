rm(list = ls())
library(INLA)
#set.seed(1)
inla.setOption(pardiso.license="~/.pardiso.lic")
n_person = 2000
n_each_person = 20
nL = n_each_person*n_person
age = sample(x = 15:75,size = n_person,replace = T)
age = rnorm(n_person)
age = scale(age)
shape = 1
ui = rnorm(n_person)
ui = scale(ui)
#ui = 0
eta_event1 = .1 + .1*age + .1*ui
eta_event2 = .2 + .2*age + .2*ui
eta_event3 = .3 + .3*age + .3*ui
rate1 = exp(eta_event1)
rate2 = exp(eta_event2)
rate3 = exp(eta_event3)
T1 = rexp(n = n_person,rate = rate1)
T2 = rexp(n = n_person,rate = rate2)
T3 = rexp(n = n_person,rate = rate3)
Time = numeric(n_person)
Event = numeric(n_person)
Event1 = rep(0,n_person)
Event2 = rep(0,n_person)
Event3 = rep(0,n_person)
time_longitudinal = c()
for(i in 1:n_person){
    Event[i] = sample(x = 1:3,size = 1)
    if(Event[i] == 1){
        Time[i] = T1[i]
        Event1[i] = 1
    }else if(Event[i] == 2){
        Time[i] = T2[i]
        Event2[i] = 1
    }else if(Event[i] == 3){
        Time[i] = T3[i]
        Event3[i] = 1
    }
    time_longitudinal = c(time_longitudinal,sort(runif(n = n_each_person,min = 0,max = Time[i])))
}

eta_longitudinal = -0.4 - .4*rep(age,each = n_each_person)  - .4*time_longitudinal + rep(ui,each = n_each_person)
lambda_longitudinal = exp(eta_longitudinal)
y_longitudinal = rpois(n = nL,lambda = lambda_longitudinal)

data = list(
    y_longitudinal = c(y_longitudinal,rep(NA,3*n_person))
)






