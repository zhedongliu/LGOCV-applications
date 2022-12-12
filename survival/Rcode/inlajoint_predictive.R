rm(list = ls())
library(INLAjoint)
library(pracma)
source("~/r-inla/rinla/R/likelihood.R")
inla.setOption(pardiso.license="~/.pardiso.lic")
set.seed(10)
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
time_longitudinal1 = c()
time_longitudinal2 = c()
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
    time_longitudinal1 = c(time_longitudinal1,sort(runif(n = n_each_person,min = 0,max = Time[i])))
    time_longitudinal2 = c(time_longitudinal2,sort(runif(n = n_each_person,min = 0,max = Time[i])))
}

eta_longitudinal1 = .4 + .4*rep(age,each = n_each_person)  + .4*time_longitudinal1 + rep(ui,each = n_each_person)
eta_longitudinal2 = .5 + .5*rep(age,each = n_each_person)  + .5*time_longitudinal2 + .5*rep(ui,each = n_each_person)

lambda_longitudinal1 = exp(eta_longitudinal1)
y_longitudinal1 = rpois(n = nL,lambda = lambda_longitudinal1)

lambda_longitudinal2 = exp(eta_longitudinal2)
y_longitudinal2 = rpois(n = nL,lambda = lambda_longitudinal2)

Event1[1] = 0
Event2[1] = 0
Event3[1] = 0

y_surv1 = inla.surv(time = Time,event = Event1)
y_surv2 = inla.surv(time = Time,event = Event2)
y_surv3 = inla.surv(time = Time,event = Event3)




LongDat1 = data.frame("heart" = y_longitudinal1,"time" = time_longitudinal1,
                      "age" = rep(age,each = n_each_person), "idd" = rep(1:n_person,each = n_each_person))

LongDat2 = data.frame("blood" = y_longitudinal2,"time" = time_longitudinal2,
                       "age2" = rep(age,each = n_each_person), "idd" = rep(1:n_person,each = n_each_person))

SurvDat = data.frame("age" = age, "idd" = 1:n_person)

res = joint(formLong = list(heart ~ time + age + (1|idd),
                             blood ~ time + age2 + (1|idd)),
             formSurv = list(y_surv1 ~ age,
                             y_surv2 ~ age,
                             y_surv3 ~ age),
             id="idd", 
             timeVar= "time", 
             family = c("poisson","poisson"),
             basRisk = c("weibullsurv", "weibullsurv", "weibullsurv"),
             dataLong = list(LongDat1,LongDat2),
             dataSurv = SurvDat, 
             assoc= list(c("SRE_ind","SRE_ind","SRE_ind"),c("SRE_ind","SRE_ind","SRE_ind")),
             control=list(likelihood.info = TRUE))


#result              an INLAjoint object.
#reference_time      the time point at which you make the prediction. Default is 0.
#end_time            the end point of your prediction. NULL means you don't have a end time point. Default is NULL.
#scores              desired score you want compute. "all" means logarithmic score, Brier score and ROC-AUC will be computed. Or any one of them. "Log" is logarithmic score, "Brier" is Brier score and "AUC" is ROC-AUC. Default is "all".
#l_ROC               desired length of the sequence of criteria for ROC. Default is 30.
#sample_size         desired sample_size of Monte Carlo integration. 1 means "empirical bayes", where mode will be used.

INLAjoint_prediction_scores = function(result = res, reference_time = 0, prediction_time = NULL, scores = "all", l_ROC = 30, sample_size = 1){
    event_num = length(result$basRisk)
    if(event_num == 0) stop("No suvival object exists.")
    n_data = length(result$.args$offset)
    n_individual = sum(!is.na(unique(result$.args$data[[result$id]])))
    survidx = lapply(1:n_individual, function(x) which(result$.args$data[[result$id]]==x))
    longi_num = length(result$famLongi)
    key_name = character(longi_num)
    time_name = character(longi_num)
    #is this surve1time fiexed?
    time = result$.args$data$surv1time[!is.na(result$.args$data$surv1time)]
    
    for(longi_type in 1:longi_num){
        pattern = paste0("^ID.*L",longi_type,"$")
        key_name[longi_type] = grep(pattern = pattern,x = names(result$.args$data),value = TRUE)
        time_name[longi_type] = paste0(result$timeVar,"_L",longi_type)
    }
    
    longitudinal_offset = length(result$.args$family) - event_num
    groups = vector(mode = "list",length = n_data)
    event = numeric(n_individual)
    for(individual_idx in 1:n_individual){
        longi_idx = c()
        for(longi_type in 1:longi_num){
            longi_idx_temp = which(result$.args$data[[key_name[longi_type]]] == individual_idx)
            longi_idx_temp = longi_idx_temp[which(result$.args$data[[time_name[longi_type]]][longi_idx_temp] > reference_time)]
            longi_idx = c(longi_idx,longi_idx_temp)
        }
        
        for(data_idx in survidx[[individual_idx]]){
            groups[[data_idx]] = c(survidx[[individual_idx]],longi_idx)
        }
        
        for(event_idx in 1:event_num){
            if(result$.args$data$Yjoint[[longitudinal_offset + event_idx]]$time[survidx[[individual_idx]][event_idx]] != 0) event[individual_idx] = event_idx
        }
    }
    
    res_cv = inla.group.cv(result = result,groups = groups)
    # if we can sample .... now I do marginal sampling which is wrong. I simply use it to create the frame to compute the correct samples
    res_log_density = rep(NA,n_individual)
    res_prob = matrix(NA,n_individual,event_num + 1)
    res_c = seq(0,1,length.out = l_ROC)
    res_TPR = matrix(0,l_ROC,event_num+1)
    res_FPR = matrix(0,l_ROC,event_num+1)
    res_Brier = rep(0,event_num + 1)
    res_AUC = rep(0,event_num + 1)
    num_eval = 35
    args = vector(mode = "list",length = event_num)
    
    for(config_idx in 1:result$misc$configs$nconfig){
        if(result$misc$configs$config[[config_idx]]$log.posterior == 0){
            config_mode_idx = config_idx
        }
    }
    
    for(individual_idx in 1:n_individual){
        if(time[individual_idx] > reference_time){
            eta_sample = vector(mode = "list",length = event_num)
            if(sample_size == 1){
                for(event_idx in 1:event_num){
                    eta_sample[[event_idx]] = res_cv$mean[survidx[[individual_idx]][event_idx]]
                }
            }else{
                for(event_idx in 1:event_num){
                    eta_sample[[event_idx]] = rnorm(n = sample_size,mean = res_cv$mean[survidx[[individual_idx]][event_idx]],sd = res_cv$sd[survidx[[individual_idx]][event_idx]])
                }
            }
            
            for(event_idx in 1:event_num){
                args[[event_idx]] = inla.likelihood.parser(result$misc$configs$config[[config_mode_idx]]$arg.str[survidx[[individual_idx]][event_idx]])
                args[[event_idx]]$islog = TRUE
                args[[event_idx]]$linear.predictor = eta_sample[[event_idx]]
            }
            
           
            if(scores == "all" || scores == "Log"){
                log_density_i = 0
                for(event_idx in 1:event_num){
                    if(event_idx == event[individual_idx]){
                        pdf = inla.likelihood(type = "d",args = args[[event_idx]])
                        log_density_i = log_density_i + pdf(time[individual_idx])
                    }else{
                        surv_func = inla.likelihood(type = "s",args = args[[event_idx]])
                        log_density_i = log_density_i + surv_func(time[individual_idx])
                    }
                }
                res_log_density[individual_idx] = mean(log_density_i)
            }
            
            if(scores == "all" || scores == "Brier" || scores == "AUC"){
                prob_i = 0
                survive_lower = 0
                survive_upper = 0
                for(event_idx in 1:event_num){
                    surv_func = inla.likelihood(type = "s",args = args[[event_idx]])
                    survive_lower = survive_lower + surv_func(reference_time)
                    survive_upper = survive_upper + surv_func(prediction_time)
                }
                
                for(event_idx in 1:event_num){
                    args[[event_idx]]$islog = FALSE
                    cdf = inla.likelihood(type = "p",args[[event_idx]])
                    cdf_inv = inla.likelihood(type = "q",args[[event_idx]])
                    p_start = cdf(reference_time)
                    p_end = cdf(prediction_time)
                    eval_prob = t(matrix(unlist(lapply(X = 1:sample_size,FUN = function(i){seq(p_start[i],p_end[i],length.out = num_eval)})),nrow = num_eval,ncol = sample_size))
                    diff_eval = eval_prob[,2] - eval_prob[,1]
                    eval_points = cdf_inv(eval_prob)
                    func_val = 0
                    for(event_idx2 in 1:event_num){
                        if(event_idx2 != event_idx){
                            args[[event_idx]]$islog = TRUE
                            surv_func = inla.likelihood(type = "s",args = args[[event_idx2]])
                            func_val = func_val + surv_func(eval_points)
                        }
                    }
                    func_val = exp(func_val-survive_lower)
                    prob_i = rowSums(2*func_val) - func_val[,1] - func_val[,num_eval]
                    prob_i = prob_i*diff_eval*0.5
                    res_prob[individual_idx,event_idx] = mean(prob_i)
                }
                res_prob[individual_idx,event_idx + 1] = mean(exp(survive_upper - survive_lower))
                
                #res_prob[individual_idx,] = res_prob[individual_idx,]/denominator
                res_prob[individual_idx,] = res_prob[individual_idx,]/sum(res_prob[individual_idx,])
            }
        }
    }
    
    if(scores == "all" || scores == "Brier"){
        observed_unknown = time <= prediction_time & time > reference_time & event == 0
        for(event_idx in 1:event_num){
            observed_positive = time <= prediction_time & time > reference_time & event == event_idx
            observed_positive[observed_unknown] = NA
            res_Brier[event_idx] = mean((observed_positive - res_prob[,event_idx])^2,na.rm = TRUE)
        }
        observed_positive = time > prediction_time
        observed_positive[observed_unknown] = NA
        res_Brier[event_num + 1] = mean((observed_positive - res_prob[,event_num + 1])^2,na.rm = TRUE)
    }
    
    if(scores == "all" || scores == "AUC"){
        observed_unknown = time <= prediction_time & time > reference_time & event == 0
        for(event_idx in 1:event_num){
            observed_positive = time <= prediction_time & time > reference_time & event == event_idx
            observed_positive[observed_unknown] = NA
            observed_negative = time > prediction_time | (time > reference_time & event != event_idx)
            observed_negative[observed_unknown] = NA
            positive_fold = which(observed_positive)
            negative_fold = which(observed_negative)
            for(individual_idx in positive_fold){
                res_AUC[event_idx] = res_AUC[event_idx] + sum(res_prob[individual_idx,event_idx] > res_prob[negative_fold,event_idx])
                for(i in 1:l_ROC){
                    model_decision = res_prob[,event_idx] > res_c[i]
                    res_TPR[i,event_idx] = sum(model_decision & observed_positive,na.rm = TRUE)/sum(observed_positive,na.rm = TRUE)
                    res_FPR[i,event_idx] = sum(model_decision & observed_negative,na.rm = TRUE)/sum(observed_negative,na.rm = TRUE)
                }
            }
            res_AUC[event_idx] = res_AUC[event_idx]/length(positive_fold)/length(negative_fold)
        }
        
        observed_positive = time > prediction_time
        observed_positive[observed_unknown] = NA
        observed_negative = time <= prediction_time & time > reference_time
        observed_negative[observed_unknown] = NA
        positive_fold = which(observed_positive)
        negative_fold = which(observed_negative)
        for(individual_idx in positive_fold){
            res_AUC[event_num+1] = res_AUC[event_num+1] + sum(res_prob[individual_idx,event_num+1] > res_prob[negative_fold],na.rm = TRUE)
            for(i in 1:l_ROC){
                model_decision = res_prob[,event_num+1] > res_c[i]
                res_TPR[i,event_num+1] = sum(model_decision & observed_positive,na.rm = TRUE)/sum(observed_positive,na.rm = TRUE)
                res_FPR[i,event_num+1] = sum(model_decision & observed_negative,na.rm = TRUE)/sum(observed_negative,na.rm = TRUE)
            }
        }
        res_AUC[event_num+1] = res_AUC[event_num+1]/length(positive_fold)/length(negative_fold)
    }
    
    return(list(log_density = res_log_density, Brier_score = res_Brier, AUC = res_AUC, prob = res_prob, ROC = list(threshold = res_c, TPR = res_TPR, FPR = res_FPR)))
}

plot_ROC = function(pred_res){
    event_num = dim(pred_res$ROC$TPR)[2] - 1
    for(event_idx in 1:event_num){
        plot(pred_res$ROC$FPR[,event_idx],pred_res$ROC$TPR[,event_idx],type="b",main = paste0("ROC for event ", event_idx),xlab = "FPR",ylab = "TPR",ylim = c(0,1),xlim=c(0,1))
        abline(a = 0, b =1)
    }
    
    plot(pred_res$ROC$FPR[,event_num+1],pred_res$ROC$TPR[,event_num+1],type="b",main = paste0("ROC for no event happen"),xlab = "FPR",ylab = "TPR",ylim = c(0,1),xlim=c(0,1))
    abline(a = 0, b =1)
}

pred_res = INLAjoint_prediction_scores(result = res,reference_time = .1,prediction_time = .5,scores = "all",l_ROC = 100,sample_size = 100)
plot_ROC(pred_res)
