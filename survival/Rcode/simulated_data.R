library(INLA)
inla.setOption(pardiso.license="~/.pardiso.lic")
num_person = 1000
shape = 3
data_survival = data.frame(id = 1:num_person,time = NA,event = 1,scale = NA,personal_intercept = NA, personal_slope = NA)
censor_time = 5
rho = 0.5
for (i in 1:num_person){
    data_survival$personal_intercept[i] = rnorm(1,sd = .1)
    data_survival$personal_slope[i] = rho*data_survival$personal_intercept[i] + rnorm(1,sd = .1)
    data_survival$scale[i] = exp(1.5+3*data_survival$personal_intercept[i]+5*data_survival$personal_slope[i])
    data_survival$time[i] = rweibull(1,scale = data_survival$scale[i],shape = shape)
    if(data_survival$time[i] >= censor_time){
        data_survival$time[i] = censor_time + rnorm(1,sd = .01)
        data_survival$event[i] = 0
    }
}
time_possible = seq(0,max(floor(data_survival$time)),1)
intercept = 1
trend = 0.1

id_longi = c()
time_longi = c()
y_longi = c()
for(i in 1:num_person){
    time_i = time_possible[time_possible<=data_survival$time[i]]
    n_i = length(time_i)
    y_i = intercept + trend*time_i + data_survival$personal_intercept[i] + data_survival$personal_slope[i]*time_i + rnorm(n_i,sd = .1)
    y_longi = c(y_longi,y_i)
    time_longi = c(time_longi,time_i)
    id_longi = c(id_longi,rep(i,n_i))
}
data_longi = data.frame(id = id_longi,time = time_longi,y = y_longi)
NL = length(data_longi$id)
NS = length(data_survival$id)



# plot(x=NA,y=NA,ylim = c(min(data_longi$y),max(data_longi$y)),xlim=c(0,max(data_survival$time)),xlab="time",ylab="readings")
# for (i in 1:num_person){
#     lines(data_longi$time[data_longi$id==i],data_longi$y[data_longi$id==i],col=i)
#     if(data_survival$event[i] == 1){points(data_survival$time[i],data_longi$y[max(which(data_longi$id==i))])}else{
#         points(data_survival$time[i],data_longi$y[max(which(data_longi$id==i))],pch = 4)
#     }
#     
#     lines(c(data_longi$time[max(which(data_longi$id==i))],data_survival$time[i]),c(data_longi$y[max(which(data_longi$id==i))],data_longi$y[max(which(data_longi$id==i))]),lty = 2)
# }

plot(x=NA,y=NA,ylim = c(min(data_longi$y),max(data_longi$y)),xlim=c(0,max(data_survival$time)),xlab="time",ylab="readings")
for (i in sample(num_person,size = min(10,num_person))){
    lines(data_longi$time[data_longi$id==i],data_longi$y[data_longi$id==i],col=i)
    if(data_survival$event[i] == 1){points(data_survival$time[i],data_longi$y[max(which(data_longi$id==i))])}else{
        points(data_survival$time[i],data_longi$y[max(which(data_longi$id==i))],pch = 4)
    }
    
    lines(c(data_longi$time[max(which(data_longi$id==i))],data_survival$time[i]),c(data_longi$y[max(which(data_longi$id==i))],data_longi$y[max(which(data_longi$id==i))]),lty = 2)
}

survout = inla.surv(time = c(rep(NA,NL),data_survival$time),event = c(rep(NA,NL),data_survival$event))

dataINLA = list(total_interceptlongi = c(rep(1,NL),rep(NA,NS)),
                      personal_intercetlongi = c(data_longi$id,rep(NA,NS)),
                      personal_slopelongi = c(data_longi$id+NS,rep(NA,NS)),
                      time = c(data_longi$time,rep(NA,NS)),
                      total_interceptsurv = c(rep(NA,NL),rep(1,NS)),
                      personal_interceptsurv = c(rep(NA,NL),data_survival$id),
                      personal_slopesurv = c(rep(NA,NL),data_survival$id+NS),
                      Yjoint = list(c(data_longi$y,rep(NA,NS)),survout))

formula = Yjoint ~ -1+ total_interceptlongi + time + f(personal_intercetlongi,model = "iidkd",order=2,n = NS*2) + f(personal_slopelongi,time,copy = "personal_intercetlongi") +
                   total_interceptsurv + f(personal_interceptsurv,copy = "personal_intercetlongi",hyper = list(beta=list(fixed = FALSE))) + f(personal_slopesurv,copy = "personal_intercetlongi",hyper =  list(beta=list(fixed = FALSE)))

res = inla(formula = formula, 
               data = dataINLA, 
               family=c("gaussian", "weibullsurv"),
               inla.mode = "experimental",
               control.family=list(list(), list(variant=1)),
               control.inla = list(int.strategy="eb"),
               control.compute = list(config = T,smtp="pardiso"),verbose = T)



lgocv = inla.group.cv(res,num.level.sets = 3,verbose = F)

plot(lgocv$cv)
points(lgocv$cv,col="red")
points(lgocv$cv,col="blue")

n = length(lgocv$mean)
node = c()
groups = list()
for (i in 1:n){
    node = c(node,lgocv$groups[[i]]$idx)
}
length(unique(node))
m = matrix(NA,length(lgocv$mean),length(lgocv$mean))
for(i in 1:n){
    l = length(lgocv$groups[[i]]$idx)
    for (ii in 1:l){
        for(jj in 1:l){
            m[lgocv$groups[[i]]$idx[ii],lgocv$groups[[i]]$idx[jj]] = 1
        }
    }
}
(sum(!is.na(m))-sum(!is.na(diag(m))))/2
m[is.na(m)] = 0
image(m)
# res$summary.fixed
# res$summary.hyperpar
# 
# 
# 
# 
# 
# 
# MC_samples <- inla.iidkd.sample(10^4, res, "personal_intercetlongi", return.cov=TRUE) # for random-effects covariance terms
# VarCov <- matrix(unlist(MC_samples), nrow = 2^2)
# VarCovMeans <- matrix(rowMeans(VarCov),2,2)
# VarCovSD <- matrix(apply(VarCov,1,sd),2,2)
# VarCov025 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.025)),2,2)
# VarCov975 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.975)),2,2)
# 
# 
# sqrt((VarCovMeans[2,2] - 0.1^2)/VarCovMeans[1,1])
# VarCovMeans
# VarCovMeans/outer(sqrt(diag(VarCovMeans)),sqrt(diag(VarCovMeans)))
# 
# a = exp(res$summary.hyperpar$mode[3])
# b = exp(res$summary.hyperpar$mode[4])
# p = res$summary.hyperpar$mode[5]
# 2*exp(p)/(1+exp(p)) - 1
