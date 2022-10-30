# reference with frailtypack
library(frailtypack)
data(colorectal)
data(colorectalLongi)

# Survival data preparation - only terminal events
colorectalSurv <- subset(colorectal, new.lesions == 0)
# Box-cox back transformation (lambda=0.3) and apply logarithm (with a 1 unit shift)
colorectalLongi$Yo <- (colorectalLongi$tumor.size*0.3+1)^(1/0.3)
colorectalLongi$Y <- log(colorectalLongi$Yo+1) # log transformation with shift=1
# modLongi <- longiPenal(Surv(time1, state) ~ treatment + prev.resection,
#                        Y ~  year*treatment + prev.resection, colorectalSurv, data.Longi =colorectalLongi,
#                        random = c("1", "year"),  id = "id", link = "Random-effects",
#                        hazard = "Weibull", method.GH = "Pseudo-adaptive", n.nodes = 32)
# saveRDS(object = modLongi,file = "~/Dropbox/phdproject/survival/modLongi.RDS")
modLongi = readRDS("~/Dropbox/phdproject/survival/modLongi.RDS")
# computation on the same dataset
time <- .0001
epoce <- epoce(modLongi,time)

print(epoce)
plot(epoce, type = "cvpol")


# same model with INLA
library(INLA)
NL <- nrow(colorectalLongi) # number of repeated measurement (longitudinal continuous outcome)
NS <- nrow(colorectalSurv) # number of survival times (i.e., number of individuals)

# survival outcome
SurvOut <- inla.surv(time = c(rep(NA, NL), colorectalSurv$time1), event = c(rep(NA, NL), colorectalSurv$state))
# create a dataset with two likelihoods for INLA
dataINLA = list(
  IDrandominterceptLongi = c(colorectalLongi$id, rep(NA, NS)),
  IDrandomslopeLongi = c(NS+colorectalLongi$id, rep(NA, NS)),
  InterceptLongi = c(rep(1, NL), rep(NA, NS)),
  YearLongi = c(colorectalLongi$year, rep(NA, NS)),
  TrtLongi = c(colorectalLongi$treatment, rep(NA, NS)),
  Slope.X.TrtLongi = c(colorectalLongi$year*(as.integer(colorectalLongi$treatment)-1), rep(NA, NS)),
  PrevresLongi = c(colorectalLongi$prev.resection, rep(NA, NS)),
  IDrandominterceptSurv = c(rep(NA, NL), colorectalSurv$id),
  IDrandomslopeSurv = c(rep(NA, NL), NS+colorectalSurv$id),
  InterceptSurv = c(rep(NA, NL), rep(1, NS)),
  TrtSurv = c(rep(NA, NL), colorectalSurv$treatment),
  PrevresSurv = c(rep(NA, NL), colorectalSurv$prev.resection),
  Yjoint = list(c(colorectalLongi$Y, rep(NA, NS)), SurvOut)
)

formJoint = Yjoint ~ -1 + InterceptLongi + YearLongi + TrtLongi + Slope.X.TrtLongi + PrevresLongi + #longitudinal part of the formula (fixed effects)
                          InterceptSurv + TrtSurv + PrevresSurv + # survival part of the formula (fixed effects)
                            f(IDrandominterceptLongi, model = "iidkd", order=2, n=NS*2) + # random intercept
                            f(IDrandomslopeLongi, YearLongi, copy="IDrandominterceptLongi") + # random slope
                            f(IDrandominterceptSurv, copy="IDrandominterceptLongi", fixed=FALSE) + # random intercept shared in survival model (scaled by beta)
                            f(IDrandomslopeSurv, copy="IDrandominterceptLongi", fixed=FALSE) # random slope shared in survival model (scaled by beta)
modINLA = inla(formula = formJoint, 
               data=dataINLA,
               inla.mode = "experimental",
               family=c("gaussian", "weibullsurv"),
               control.family=list(list(), list(variant=1)), 
               control.inla=list(int.strategy="eb"),
               control.compute = list(config = T,control.gcpo = list(enable = TRUE)))

modINLA$summary.fixed
modINLA$summary.hyperpar


MC_samples <- inla.iidkd.sample(10^4, modINLA, "IDrandominterceptLongi", return.cov=TRUE) # for random-effects covariance terms
VarCov <- matrix(unlist(MC_samples), nrow = 2^2)
VarCovMeans <- matrix(rowMeans(VarCov),2,2)
VarCovSD <- matrix(apply(VarCov,1,sd),2,2)
VarCov025 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.025)),2,2)
VarCov975 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.975)),2,2)


loocv = inla.group.cv(modINLA)
plot(epoce$IndivContrib[,1][colorectalSurv$state==1],log(loocv$cv[which(!is.na(SurvOut$time))][colorectalSurv$state==1]),xlab = "EPOCE(time = 0.0001)",ylab = "LOOCV_LOG")
abline(a=0,b=1)

plot(epoce$IndivContrib[,1][colorectalSurv$state==0],log(loocv$cv[which(!is.na(SurvOut$time))][colorectalSurv$state==0]),xlab = "EPOCE(time = 0.0001)",ylab = "LOOCV_LOG")
abline(a=0,b=1)




alpha = unname(exp(0.1*modINLA$misc$configs$config[[1]]$theta[2]))
eta_mean_all = modINLA$gcpo$mean
eta_sd_all = modINLA$gcpo$sd
eta_mean = eta_mean_all[which(!is.na(SurvOut$time))][colorectalSurv$state==0]
eta_sd = eta_sd_all[which(!is.na(SurvOut$time))][colorectalSurv$state==0]
loocv_censored = numeric(length(modINLA$summary.linear.predictor$mean[which(!is.na(SurvOut$time))][colorectalSurv$state==0]))


for (point.interest in 1:29){
    y.here = SurvOut$lower[which(!is.na(SurvOut$time))][colorectalSurv$state==0][point.interest]
    xx = seq(eta_mean[1] - 20*eta_sd[point.interest],eta_mean[point.interest] + 20*eta_sd[point.interest],length.out = 1000)
    #yy = pweibull(y.here,shape = alpha,scale = (exp(xx))^(-1/alpha),lower.tail = FALSE) * dnorm(xx,mean = eta_mean[point.interest],sd = eta_sd[1])
    yy = dweibull(y.here,shape = alpha,scale = (exp(xx))^(-1/alpha)) * dnorm(xx,mean = eta_mean[point.interest],sd = eta_sd[1])
    plot(xx,yy)
    loocv_censored[point.interest] = pracma::trapz(xx,yy)
}

plot(log(loocv_censored),epoce$IndivContrib[,1][colorectalSurv$state==0])
plot(log(loocv$cv[which(!is.na(SurvOut$time))][colorectalSurv$state==0]),epoce$IndivContrib[,1][colorectalSurv$state==0])
plot(log(loocv_censored),log(loocv$cv[which(!is.na(SurvOut$time))][colorectalSurv$state==0]))
abline(a=0,b=1)
