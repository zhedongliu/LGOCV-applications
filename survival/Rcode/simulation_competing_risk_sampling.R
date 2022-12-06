rm(list = ls())
library(INLA)
library(pracma)
source("~/r-inla/rinla/R/likelihood.R")
inla.setOption(pardiso.license="~/.pardiso.lic")
n_person = 1000
n_each_person = 5
nL = n_each_person*n_person
age = sample(x = 15:75,size = n_person,replace = T)
age = scale(age)
shape = 3
theta = 10*log(shape)
ui = rnorm(n_person)
ui = scale(ui)

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
    Event[i] = which.min(c(T1[i],T2[i],T3[i]))
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

eta_longitudinal = .4 + .4*rep(age,each = n_each_person)  + .4*time_longitudinal + rep(ui,each = n_each_person)
lambda_longitudinal = exp(eta_longitudinal)
y_longitudinal = rpois(n = nL,lambda = lambda_longitudinal)

y_longitudinal = c(y_longitudinal,rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
y_surv1 = inla.surv(time = c(rep(NA,nL),Time,rep(NA,n_person),rep(NA,n_person)),event = c(rep(NA,nL),Event1,rep(NA,n_person),rep(NA,n_person)))
y_surv2 = inla.surv(time = c(rep(NA,nL),rep(NA,n_person),Time,rep(NA,n_person)),event = c(rep(NA,nL),rep(NA,n_person),Event2,rep(NA,n_person)))
y_surv3 = inla.surv(time = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),Time),event = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),Event3))

id_longitudinal = c(rep(1:n_person,each = n_each_person),rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
id_surv1 = c(rep(NA,nL),1:n_person,rep(NA,n_person),rep(NA,n_person))
id_surv2 = c(rep(NA,nL),rep(NA,n_person),1:n_person,rep(NA,n_person))
id_surv3 = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),1:n_person)

age_longitudinal = c(rep(age,each = n_each_person),rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
age_surv1 = c(rep(NA,nL),age,rep(NA,n_person),rep(NA,n_person))
age_surv2 = c(rep(NA,nL),rep(NA,n_person),age,rep(NA,n_person))
age_surv3 = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),age)

mu = as.factor(c(rep(1,nL),rep(2,n_person),rep(3,n_person),rep(4,n_person)))

time_longitudinal_var = c(time_longitudinal,rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))

data_list = list(y = list(y_longitudinal,y_surv1,y_surv2,y_surv3),
                 id_longitudinal = id_longitudinal,
                 mu = mu,
                 time_longitudinal = time_longitudinal_var,
                 age_longitudinal = age_longitudinal,
                 id_surv1 = id_surv1,
                 id_surv2 = id_surv2,
                 id_surv3 = id_surv3,
                 age_surv1 = age_surv1,
                 age_surv2 = age_surv2,
                 age_surv3 = age_surv3)


hyper_fixed = FALSE
hyper_initial = c(theta,theta,theta,0,.1,.2,.3)

formula = y ~ -1  + mu +
    time_longitudinal +
    age_longitudinal + age_surv1 + age_surv2 + age_surv3 +
    f(id_longitudinal,model = "iid",hyper = list(theta = list(initial = hyper_initial[4],fixed = hyper_fixed))) +
    f(id_surv1,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[5]))) +
    f(id_surv2,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[6]))) +
    f(id_surv3,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[7])))

res = inla(formula, 
           family = c("poisson","weibullsurv","weibullsurv","weibullsurv"),
           data = data_list, 
           control.family = list(list(),list(variant = 1,initial = hyper_initial[1], fixed = hyper_fixed),list(variant = 1,initial = hyper_initial[2], fixed = hyper_fixed),list(variant = 1,initial = hyper_initial[3], fixed = hyper_fixed)),
           inla.mode = "experimental",
           verbose = FALSE,
           control.compute = list(smtp="pardiso",config = TRUE,likelihood.info = TRUE),
           safe = TRUE)


posterior_sampling = inla.posterior.sample(n = 1000,result = res)









