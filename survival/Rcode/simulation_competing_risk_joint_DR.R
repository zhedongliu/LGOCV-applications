#rm(list = ls())
library(INLAjoint)
library(pracma)
#source("~/r-inla/rinla/R/likelihood.R")
inla.setOption(pardiso.license="~/.pardiso.lic")
n_person = 10
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
#
# y_longitudinal = c(y_longitudinal,rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
# y_surv1 = inla.surv(time = c(rep(NA,nL),Time,rep(NA,n_person),rep(NA,n_person)),event = c(rep(NA,nL),Event1,rep(NA,n_person),rep(NA,n_person)))
# y_surv2 = inla.surv(time = c(rep(NA,nL),rep(NA,n_person),Time,rep(NA,n_person)),event = c(rep(NA,nL),rep(NA,n_person),Event2,rep(NA,n_person)))
# y_surv3 = inla.surv(time = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),Time),event = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),Event3))
#
# id_longitudinal = c(rep(1:n_person,each = n_each_person),rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
# id_surv1 = c(rep(NA,nL),1:n_person,rep(NA,n_person),rep(NA,n_person))
# id_surv2 = c(rep(NA,nL),rep(NA,n_person),1:n_person,rep(NA,n_person))
# id_surv3 = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),1:n_person)
#
# age_longitudinal = c(rep(age,each = n_each_person),rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
# age_surv1 = c(rep(NA,nL),age,rep(NA,n_person),rep(NA,n_person))
# age_surv2 = c(rep(NA,nL),rep(NA,n_person),age,rep(NA,n_person))
# age_surv3 = c(rep(NA,nL),rep(NA,n_person),rep(NA,n_person),age)
#
# mu = as.factor(c(rep(1,nL),rep(2,n_person),rep(3,n_person),rep(4,n_person)))
#
# time_longitudinal_var = c(time_longitudinal,rep(NA,n_person),rep(NA,n_person),rep(NA,n_person))
#
# data_list = list(y = list(y_longitudinal,y_surv1,y_surv2,y_surv3),
#                  id_longitudinal = id_longitudinal,
#                  mu = mu,
#                  time_longitudinal = time_longitudinal_var,
#                  age_longitudinal = age_longitudinal,
#                  id_surv1 = id_surv1,
#                  id_surv2 = id_surv2,
#                  id_surv3 = id_surv3,
#                  age_surv1 = age_surv1,
#                  age_surv2 = age_surv2,
#                  age_surv3 = age_surv3)
#
#
# hyper_fixed = FALSE
# hyper_initial = c(theta,theta,theta,0,.1,.2,.3)
#
# formula = y ~ -1  + mu +
#     time_longitudinal +
#     age_longitudinal + age_surv1 + age_surv2 + age_surv3 +
#     f(id_longitudinal,model = "iid",hyper = list(theta = list(initial = hyper_initial[4],fixed = hyper_fixed))) +
#     f(id_surv1,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[5]))) +
#     f(id_surv2,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[6]))) +
#     f(id_surv3,copy = "id_longitudinal",hyper = list(theta = list(initial = hyper_initial[7])))
#
# res = inla(formula,
#            family = c("poisson","weibullsurv","weibullsurv","weibullsurv"),
#            data = data_list,
#            control.family = list(list(),list(variant = 1,initial = hyper_initial[1], fixed = hyper_fixed),list(variant = 1,initial = hyper_initial[2], fixed = hyper_fixed),list(variant = 1,initial = hyper_initial[3], fixed = hyper_fixed)),
#            inla.mode = "experimental",
#            verbose = FALSE,
#            control.compute = list(smtp="pardiso",config = TRUE,likelihood.info = TRUE),
#            safe = TRUE)

y_surv1 = inla.surv(time = Time,event = Event1)
y_surv2 = inla.surv(time = Time,event = Event2)
y_surv3 = inla.surv(time = Time,event = Event3)

LongDat <- data.frame("y_longitudinal" = y_longitudinal,"time_longitudinal" = time_longitudinal,
                      "age" = rep(age,each = n_each_person), "ID" = rep(1:n_person,each = n_each_person))
SurvDat <- data.frame("age" = age, "ID" = 1:n_person)
res <- joint(formLong = y_longitudinal ~ time_longitudinal + age + (1|ID), data=LongDat,
              formSurv = list(y_surv1 ~ age,
                              y_surv2 ~ age,
                              y_surv3 ~ age),
              id="ID", timeVar="time_longitudinal", family="poisson",
              basRisk=c("weibullsurv", "weibullsurv", "weibullsurv"),
              dataSurv = SurvDat, assoc=c("SRE_ind", "SRE_ind", "SRE_ind"),
             control=list(likelihood.info = TRUE))

friends <- lapply(res$.args$data[[res$id]], function(x) which(res$.args$data[[res$id]]==x))
individual_data <- lapply(1:n_person, function(x) which(res$.args$data[[res$id]]==x))

res$timeVar

args = list()
event_num = length(grep("surv", res$.args$family))
config_idx = 1
individual_idx = 1 # input of the function
for(event_idx in 1:event_num){
  args[[event_idx]] = inla.likelihood.parser(res$misc$configs$config[[config_idx]]$arg.str[which(res$.args$data[[res$id]]==individual_idx)[event_idx]])
}

sample_i = function(n,args){
  time_event = vector(mode = "list",length = event_num)
  for(event_idx in 1:event_num){
    sampler = inla.likelihood(type = "r",args = args[[event_idx]])
    time_event[[event_idx]] = sampler(n)
  }
  Time = numeric(n)
  Event = numeric(n)
  for(sample_idx in 1:n){
    min_time = Inf
    for(event_idx in 1:event_num){
      if(time_event[[event_idx]][sample_idx] < min_time){
        Event[sample_idx] = event_idx
        min_time = time_event[[event_idx]][sample_idx]
      }
      Time[sample_idx] = min_time
    }
  }

  return(list(Time = Time,Event = Event))
}
sample_i(1,args)


#pi(time = t, event = j)
density_i = function(time,event,args,islog = FALSE){
  res = numeric(length(density_i))
  for(event_idx in 1:event_num){
    if(event_idx == event){
      args[[event_idx]]$islog = TRUE
      pdf = inla.likelihood(type = "d",args = args[[event_idx]])
      res = res + pdf(time)
    }else{
      args[[event_idx]]$islog = TRUE
      surv_func = inla.likelihood(type = "s",args = args[[event_idx]])
      res = res + surv_func(time)
    }
  }
  if(!islog){res = exp(res)}
  return(res)
}


#P(s< time <= t, event = j)
prob_i = function(lower,upper,event,args,islog = FALSE){
  num_eval = 35
  cdf = inla.likelihood(type = "p",args[[event]])
  eval_prob = seq(cdf(lower),cdf(upper),length.out = num_eval)
  eval_points = inla.likelihood(type = "q",args[[event]])(eval_prob)
  func_val = numeric(num_eval)
  for(event_idx in 1:event_num){
    if(event_idx != event){
      args[[event_idx]]$islog = TRUE
      surv_func = inla.likelihood(type = "s",args = args[[event_idx]])
      func_val = func_val + surv_func(eval_points)
    }
  }
  func_val = exp(func_val)
  res = sum(2*func_val) - func_val[1] - func_val[num_eval]
  res = res * (eval_prob[2]-eval_prob[1]) * 0.5
  if(islog){res = log(res)}
  return(res)
}
#P(time>=s, event = j)
surv_i = function(lower,args,islog = FALSE){
  res = 0
  for(event_idx in 1:event_num){
    args[[event_idx]]$islog = TRUE
    surv_func = inla.likelihood(type = "s",args = args[[event_idx]])
    res = res + surv_func(lower)
  }
  if(!islog){res = exp(res)}
  return(res)
}

conditional_prob_i = function(lower,upper,args,islog = FALSE){
  p = numeric(event_num+1)
  surv_lower = surv_i(lower = lower,args = args)
  for(event_idx in 1:event_num){
    p[event_idx] = prob_i(lower = lower,upper = upper,event = event_idx,args = args)
  }
  p[event_idx+1] = surv_i(lower = upper,args = args)
  p = p/surv_lower
  p = p/sum(p)
  if(islog){return(log(p))}
  return(p)
}

conditional_logarithmic_score = function(res,lower = 0){
  #get groups from somewhere
  res_cv = inla.group.cv(result = res,groups = friends)
  mode = list()
  mode$eta = res_cv$mean
  args = list()
  #get n_person from somewhere
  #get individual_data from somewhere
  #get Time from somewhere
  #get Event from somewhere
  #get mode config_idx
  config_idx = 1
  log_density = rep(NA,n_person)
  for(individual_idx in 1:n_person){
    if(Time[individual_idx] > lower & Event[individual_idx] != 0){
      for(event_idx in 1:event_num){
        args[[event_idx]] = inla.likelihood.parser(res$misc$configs$config[[config_idx]]$arg.str[individual_data[[individual_idx]][event_idx]])
        args[[event_idx]]$linear.predictor = mode$eta[individual_data[[individual_idx]][event_idx]]
      }
      log_density[individual_idx] = density_i(time = Time[individual_idx],event = Event[individual_idx],args = args,islog = TRUE) - surv_i(lower = lower,args = args,islog = TRUE)
    }
  }
  return(log_density)
}
conditional_logarithmic_score(res,0)

conditional_prob = function(res,lower,upper){
  #get groups from somewhere
  res_cv = inla.group.cv(result = res,groups = friends)
  mode = list()
  mode$eta = res_cv$mean
  args = list()
  #get n_person from somewhere
  #get individual_data from somewhere
  #get Time from somewhere
  #get Event from somewhere
  #get mode config_idx
  config_idx = 1
  prob = matrix(NA,n_person,event_num + 1)
  for(individual_idx in 1:n_person){
    for(event_idx in 1:event_num){
      args[[event_idx]] = inla.likelihood.parser(res$misc$configs$config[[config_idx]]$arg.str[individual_data[[individual_idx]][event_idx]])
      args[[event_idx]]$linear.predictor = mode$eta[individual_data[[individual_idx]][event_idx]]
    }
    prob[individual_idx,] = conditional_prob_i(lower = lower,upper = upper, args = args)
  }
  return(prob)
}

prob = conditional_prob(res = res,0,0.5)

AUC_compute = function(res = res,lower,upper,prob = NULL){
  if(is.null(prob)){prob = conditional_prob(res = res,lower,upper)}
  AUC = numeric(event_num+1)
  for(event_idx in 1:event_num){
    true_fold = which(Time <= upper & Time > lower & Event == event_idx)
    false_fold = which(Time > upper | (Time > lower & Event != event_idx))
    for(individual_idx in true_fold){
      AUC[event_idx] = AUC[event_idx] + sum(prob[individual_idx,event_idx] > prob[false_fold,event_idx])
    }
    AUC[event_idx] = AUC[event_idx]/length(true_fold)/length(false_fold)
  }
  true_fold = which(Time > upper)
  false_fold = which(Time <= upper & Time > lower)
  for(individual_idx in true_fold){
    AUC[event_num+1] = AUC[event_num+1] + sum(prob[individual_idx,event_num+1] > prob[false_fold,event_num+1])
  }
  AUC[event_num+1] = AUC[event_num+1]/length(true_fold)/length(false_fold)
  return(AUC)
}

ROC_compute = function(res,lower,upper,n = 100,prob = NULL){
  if(is.null(prob)){prob = conditional_prob(res = res,lower,upper)}
  c = seq(1,0,length.out = n)
  TPR = matrix(NA,n,event_num+1)
  FPR = matrix(NA,n,event_num+1)
  for(event_idx in 1:event_num){
    observed_truth = Time <= upper & Time > lower & Event == event_idx
    observed_false = Time > upper | (Time > lower & Event != event_idx)
    for(i in 1:n){
      model_decision = prob[,event_idx] > c[i]
      TPR[i,event_idx] = sum(model_decision & observed_truth)/sum(observed_truth)
      FPR[i,event_idx] = sum(model_decision & observed_false)/sum(observed_false)
    }
    plot(FPR[,event_idx],TPR[,event_idx],type="b",main = paste0("ROC for event ", event_idx),xlab = "FPR",ylab = "TPR",ylim = c(0,1),xlim=c(0,1))
    abline(a = 0, b =1)
  }

  observed_truth = Time > upper
  observed_false = Time <= upper & Time > lower
  for(i in 1:n){
    model_decision = prob[,event_num+1] > c[i]
    TPR[i,event_num + 1] = sum(model_decision & observed_truth)/sum(observed_truth)
    FPR[i,event_num + 1] = sum(model_decision & observed_false)/sum(observed_false)
  }
  plot(FPR[,event_num+1],TPR[,event_num+1],type="b",main = paste0("ROC for no event happen"),xlab = "FPR",ylab = "TPR",ylim = c(0,1),xlim=c(0,1))
  abline(a = 0, b =1)
  return(list(c = c, TPR = TPR, FPR = FPR))
}


Brier_compute = function(res,lower,upper,prob = NULL){
  if(is.null(prob)){prob = conditional_prob(res = res,lower,upper)}
  Brier_sub = rep(NA,event_num + 1)
  for(event_idx in 1:event_num){
    observed_truth = Time <= upper & Time > lower & Event == event_idx
    Brier_sub[event_idx] = mean((observed_truth - prob[,event_idx])^2)
  }
  observed_truth = Time > upper
  Brier_sub[event_num + 1] = mean((observed_truth - prob[,event_num + 1])^2)
  return(Brier_sub)
}



AUC_compute(res,0,0.5,prob = prob)
ROC = ROC_compute(res,0,0.5,prob = prob)
Brier = Brier_compute(res,0,0.5,prob = prob)



