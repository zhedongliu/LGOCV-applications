library(profvis)
rm(list = ls())
library(INLA)
library(pracma)
source("~/r-inla/rinla/R/likelihood.R")
inla.setOption(pardiso.license="~/.pardiso.lic")
n_person = 100
n_each_person = 5
nL = n_each_person*n_person
age = sample(x = 15:75,size = n_person,replace = T)
age = scale(age)
shape = 3
theta = 10*log(shape)
ui = rnorm(n_person,sd = 1)
#ui = scale(ui)

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
          control.inla = list(int.strategy = "eb"),
          safe = TRUE)

s = 0 #reference time
likelihood_str = res$misc$configs$config[[1]]$arg.str
id = c(rep(1:n_person,each = n_each_person),rep(1:n_person,3))
n = length(likelihood_str)
individual_LG = list()
individual_SV = list()


for(i in 1:n_person){
    individual_LG[[i]] = numeric(0)
    individual_SV[[i]] = numeric(0)
}


for(point.interest in 1:n){
    arg = eval(parse(text = likelihood_str[point.interest]))
    if(arg$family != "inla.surv"){
        individual_LG[[id[point.interest]]] = c(individual_LG[[id[point.interest]]],point.interest)
    }else{
        individual_SV[[id[point.interest]]] = c(individual_SV[[id[point.interest]]],point.interest)
    }
}

groups = list()

for(i in 1:n){
    groups[[i]] = numeric(0)
}

for(individual.interest in 1:n_person){
    if(Time[individual.interest] > s){
        future_lg_idx = individual_LG[[individual.interest]][which(time_longitudinal[individual_LG[[individual.interest]]] > s)]
        for(point.interest in individual_SV[[individual.interest]]){
            groups[[point.interest]] = c(individual_SV[[individual.interest]],future_lg_idx)
        }
    }
}
res_cv = inla.group.cv(res = res,groups = groups)


scores_computing = function(sample){
    event_num = 3
    t = 0.5 #prediction time
    conditional_density = rep(NA,n_person)
    prob_res = matrix(NA,n_person,event_num+1)
    indicator_res = matrix(NA,n_person,event_num+1)
    brier_score_sub = matrix(NA,n_person,event_num+1)
    brier_score = rep(NA,n_person)
    args = vector(mode = "list",length = event_num)
    d_fun = vector(mode = "list",length = event_num)
    s_fun = vector(mode = "list",length = event_num)
    prob = numeric(event_num)
    for(individual.interest in 1:n_person){
        if(Time[individual.interest] >= s){
            lower = 1
            no_event = 1
            event_idx = 1
            for(nodes in individual_SV[[individual.interest]]){
                args[[event_idx]] = inla.likelihood.parser(arg_string = likelihood_str[nodes])
                args[[event_idx]]$linear.predictor = sample$eta[nodes]
                args[[event_idx]]$theta = sample$theta[event_idx]
                d_fun[[event_idx]] = inla.likelihood("d",args[[event_idx]])
                s_fun[[event_idx]] = inla.likelihood("s",args[[event_idx]])
                event_idx = event_idx + 1
            }

            for(event_idx in 1:event_num){
                lower = lower*s_fun[[event_idx]](s)
                no_event = no_event*s_fun[[event_idx]](t)
                f = function(x){
                    res = d_fun[[event_idx]](x)
                    for(i in 1:event_num){
                        if(i != event_idx){
                            res = res*s_fun[[i]](x)
                        }
                    }
                    return(res)
                }
                prob[event_idx] = quadl(f = f,xa = s,xb = t)
                if(event_idx == Event[individual.interest]) conditional_density[individual.interest] = f(Time[individual.interest])
            }
            conditional_density[individual.interest] = conditional_density[individual.interest]/lower
            prob = prob/lower
            prob_res[individual.interest,] = c(prob,no_event)
            prob_res[individual.interest,] = prob_res[individual.interest,]/sum(prob_res[individual.interest,])
            indicator = rep(NA,event_num)
            for(i in 1:event_idx){
                indicator[i] = Event[individual.interest] == i && Time[individual.interest] <= t
            }
            indicator = c(indicator,!any(indicator))
            indicator_res[individual.interest,] = indicator
            brier_score_sub[individual.interest,] = (prob_res[individual.interest,] - indicator_res[individual.interest,])^2
            brier_score[individual.interest] = sum(brier_score_sub[individual.interest,])
        }
    }

    S_AUC = rep(0,event_num + 1)
    for(major_event in 1:(event_num)){
        true_fold = which((Time <= t & Time >= s & Event == major_event))
        false_fold = which(Time >=s & (Time >t | Event != major_event ))
        p = prob_res[,major_event]
        for(i in true_fold){
            S_AUC[major_event] = S_AUC[major_event] + sum(p[i] > p[false_fold])
        }
        S_AUC[major_event] = S_AUC[major_event]/length(true_fold)/length(false_fold)
    }



    true_fold = which(Time > t)
    false_fold = which(Time >=s &Time <=t)
    p = prob_res[,event_num + 1]
    for(i in true_fold){
        S_AUC[event_num + 1] = S_AUC[event_num + 1] + sum(p[i] > p[false_fold])
    }
    S_AUC[event_num + 1] = S_AUC[event_num + 1]/length(true_fold)/length(false_fold)

    return(list(S_AUC = S_AUC,brier_score = mean(brier_score),log_conditional_density = mean(log(conditional_density)),prob_res = prob_res))
}

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

# eta_true = c(eta_longitudinal,eta_event1,eta_event2,eta_event3)
# sample = list(eta = eta_true,theta = rep(theta,3))
# scores_computing(sample)$S_AUC
sample = list(eta = res_cv$mean,theta = res$mode$theta[1:3])
score_res = scores_computing(sample)
# sample = list(eta = res$mode$x,theta = res$mode$theta[1:3])
# scores_computing(sample)$S_AUC


TPR = function(c,i,p){
    sum((p[,i]>c)*(Time <= t & Time >= s & Event == i))/sum(Time <= t & Time >= s & Event == i)
}
FPR = function(c,i,p){
    sum((p[,i]>c)*(Time > t | (Time >= s & Event != i)))/sum(Time > t | (Time >= s & Event != i))
}
c = seq(1,0,-0.001)
xx1 = unlist(lapply(c,FPR,i = 1,p = score_res$prob_res))
yy1 = unlist(lapply(c,TPR,i = 1,p = score_res$prob_res))
plot(xx1,yy1,type="l")
trapz(xx1,yy1)
score_res$S_AUC

xx2 = unlist(lapply(c,FPR,i = 1,p = p))
yy2 = unlist(lapply(c,TPR,i = 1,p = p))
plot(xx2,yy2,type="l")
trapz(xx2,yy2)
AUC

# profvis({
#     score_res = scores_computing(sample)
# })
