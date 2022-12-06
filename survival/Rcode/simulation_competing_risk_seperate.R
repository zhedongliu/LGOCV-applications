rm(list = ls())
library(INLA)
set.seed(1)
inla.setOption(pardiso.license="~/.pardiso.lic")
n_person = 10000
n_each_person = 5
nL = n_each_person*n_person
age = sample(x = 15:75,size = n_person,replace = T)
age = scale(age)
shape = 3
ui = rnorm(n_person)
ui = scale(ui)
#ui = 0
eta_event1 = .1 + .1*age + .1*ui
eta_event2 = .2 + .2*age + .2*ui
eta_event3 = .3 + .3*age + .3*ui
scale1 = exp(-eta_event1)
scale2 = exp(-eta_event2)
scale3 = exp(-eta_event3)
T1 = rweibull(n = n_person,shape = shape,scale = scale1)
T2 = rweibull(n = n_person,shape = shape,scale = scale2)
T3 = rweibull(n = n_person,shape = shape,scale = scale3)
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



id_longitudinal = rep(1:n_person,each = n_each_person)
age_longitudinal = rep(age,each = n_each_person)

eta_longitudinal = .4 - .4*age_longitudinal  - .4*time_longitudinal + rep(ui,each = n_each_person)
lambda_longitudinal = exp(eta_longitudinal)
y_longitudinal = rpois(n = nL,lambda = lambda_longitudinal)

formula_longi = y_longitudinal ~ 1 + age_longitudinal + time_longitudinal + f(id_longitudinal,model = "iid")
res_longi = inla(formula = formula_longi,
                 family = "poisson",
                 data = list(y_longitudinal = y_longitudinal,age_longitudinal = age_longitudinal,time_longitudinal = time_longitudinal,id_longitudinal = id_longitudinal))
summary(res_longi)

y_surv1 = inla.surv(time = Time,event = Event1)
formula_surv1 = y_surv1 ~ 1 + age + f(id,model = "iid")
res_surv1 = inla(formula = formula_surv1,
                 family = "weibullsurv",
                 data = list(y_surv1 = y_surv1,age = age,id = 1:n_person),
                 control.family = list(variant = 1))
summary(res_surv1)


y_surv2 = inla.surv(time = Time,event = Event2)
formula_surv2 = y_surv2 ~ 1 + age + f(id,model = "iid")
res_surv2 = inla(formula = formula_surv2,
                 family = "weibullsurv",
                 data = list(y_surv2 = y_surv2,age = age,id = 1:n_person),
                 control.family = list(variant = 1))
summary(res_surv2)

y_surv3 = inla.surv(time = Time,event = Event3)
formula_surv3 = y_surv3 ~ 1 + age + f(id,model = "iid")
res_surv3 = inla(formula = formula_surv3,
                 family = "weibullsurv",
                 data = list(y_surv3 = y_surv3,age = age,id = 1:n_person),
                 control.family = list(variant = 1))
summary(res_surv3)
