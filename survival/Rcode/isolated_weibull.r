rm(list=ls())
set.seed(1)
library(INLA)
n = 10000
shape = 3
theta = 10*log(shape)
x = rnorm(n)
scale(x)
id = 1:n


#Without iid effects
eta = 1 + 0.1*x
scale_weibull = exp(-eta)
y = rweibull(n = n,shape = shape,scale = scale_weibull)
formula = y ~ 1 + x 
res_without_iid = inla(formula = formula,
           family = "weibull",
           data = list(x=x,y=y,id = id),
           control.family = list(variant = 1,prior = "flat",param = NA),
           inla.mode = "experimental",
           verbose = F,
           safe = F)

event = (y < 10)*1L
y_surv = inla.surv(time = y,event = event)
formula_surv = y_surv ~ 1 + x
res_surv_without_iid = inla(formula = formula_surv,
                family = "weibullsurv",
                data = list(x=x,y_surv = y_surv,id = id),
                control.family = list(variant = 1),
                inla.mode = "experimental",
                verbose = F,
                safe = F)

summary(res_without_iid)
summary(res_surv_without_iid)


#With iid effects
group_size = 100
eta = 1 + 0.1*x + rep(rnorm(group_size),each = n/group_size)
scale_weibull = exp(-eta)
y = rweibull(n = n,shape = shape,scale = scale_weibull)
formula = y ~ 1 + x + f(id,model = "iid")
res = inla(formula = formula,
           family = "weibull",
           data = list(x=x,y=y,id = rep(1:group_size,each = n/group_size)),
           control.family = list(variant = 1),
           inla.mode = "experimental",
           verbose = F,
           safe = F)

event = (y < 10)*1L
y_surv = inla.surv(time = y,event = event)
formula_surv = y_surv ~ 1 + x + f(id,model = "iid")
res_surv = inla(formula = formula_surv,
                family = "weibullsurv",
                data = list(x=x,y_surv = y_surv,id = rep(1:group_size,each = n/group_size)),
                control.family = list(variant = 1),
                inla.mode = "experimental",
                verbose = F,
                safe = F)

summary(res)
summary(res_surv)








n=100
x=rnorm(n)
y=rnorm(n,x,1)
data=list(y=y,x=x)
formula=y~1+x
result=inla(formula,family="gaussian",data=data,control.family=list(prior="gaussian",param = c(0,0)),inla.mode = "experimental")
result=inla(formula,family="gaussian",data=data,control.family=list(prior="flat",param = c()),inla.mode = "experimental")


