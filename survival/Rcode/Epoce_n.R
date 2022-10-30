# reference with frailtypack
library(frailtypack)
data(colorectal)
data(colorectalLongi)

# Survival data preparation - only terminal events
colorectalSurv <- subset(colorectal, new.lesions == 0)
# Box-cox back transformation (lambda=0.3) and apply logarithm (with a 1 unit shift)
colorectalLongi$Yo <- (colorectalLongi$tumor.size*0.3+1)^(1/0.3)
colorectalLongi$Y <- log(colorectalLongi$Yo+1) # log transformation with shift=1
#colorectalLongi <- colorectalLongi[colorectalLongi$Y>0,]
modLongi <- longiPenal(Surv(time1, state) ~ treatment + prev.resection,
                       Y ~  year*treatment + prev.resection, colorectalSurv, data.Longi =colorectalLongi,
                       random = c("1", "year"),  id = "id", link = "Random-effects",
                       hazard = "Weibull", method.GH = "Pseudo-adaptive", n.nodes = 32)

# computation on the same dataset
time <- c(1, 1.5, 2, 2.5)
epoce <- epoce(modLongi,time)

print(epoce)
plot(epoce, type = "cvpol")


# same model with INLA
library(INLA)
NL <- nrow(colorectalLongi) # number of repeated measurement (longitudinal continuous outcome)
NS <- nrow(colorectalSurv) # number of survival times (i.e., number of individuals)
colorectalLongi$treatment <- as.integer(colorectalLongi$treatment)-1
colorectalLongi$prev.resection <- as.integer(colorectalLongi$prev.resection)-1
colorectalSurv$treatment <- as.integer(colorectalSurv$treatment)-1
colorectalSurv$prev.resection <- as.integer(colorectalSurv$prev.resection)-1
# survival outcome
SurvOut <- inla.surv(time = c(rep(NA, NL), colorectalSurv$time1), event = c(rep(NA, NL), colorectalSurv$state))
# create a dataset with two likelihoods for INLA
dataINLA = list(
  IDrandominterceptLongi = c(colorectalLongi$id, rep(NA, NS)),
  IDrandomslopeLongi = c(NS+colorectalLongi$id, rep(NA, NS)),
  InterceptLongi = c(rep(1, NL), rep(NA, NS)),
  SlopeLongi = c(colorectalLongi$year, rep(NA, NS)),
  TrtLongi = c(colorectalLongi$treatment, rep(NA, NS)),
  Slope.X.TrtLongi = c(colorectalLongi$year*colorectalLongi$treatment, rep(NA, NS)),
  PrevresLongi = c(colorectalLongi$prev.resection, rep(NA, NS)),
  IDrandominterceptSurv = c(rep(NA, NL), colorectalSurv$id),
  IDrandomslopeSurv = c(rep(NA, NL), NS+colorectalSurv$id),
  InterceptSurv = c(rep(NA, NL), rep(1, NS)),
  TrtSurv = c(rep(NA, NL), colorectalSurv$treatment),
  PrevresSurv = c(rep(NA, NL), colorectalSurv$prev.resection),
  Yjoint = list(c(colorectalLongi$Y, rep(NA, NS)), SurvOut)
)

formJoint = Yjoint ~ -1 + InterceptLongi + SlopeLongi + TrtLongi + Slope.X.TrtLongi + PrevresLongi + #longitudinal part of the formula (fixed effects)
                          InterceptSurv + TrtSurv + PrevresSurv + # survival part of the formula (fixed effects)
                           f(IDrandominterceptLongi, model = "iidkd", order=2, n=NS*2, constr=F, hyper=list(theta1 = list(param = c(10, 1, 1, 0)))) + # random intercept
                           f(IDrandomslopeLongi, SlopeLongi, copy="IDrandominterceptLongi") + # random slope
                           f(IDrandominterceptSurv, copy="IDrandominterceptLongi", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.01), initial = 0.1))) + # random intercept shared in survival model (scaled by beta)
                           f(IDrandomslopeSurv, copy="IDrandomslopeLongi", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.01), initial = 0.1))) # random slope shared in survival model (scaled by beta)
modINLA = inla(formula = formJoint, data=dataINLA, family=c("gaussian", "weibullsurv"),
                control.family=list(list(), list(variant=1)),
               control.inla = list(int.strategy="eb"),
               control.fixed=list(mean=0, prec=0.01, mean.intercept=0, prec.intercept=0.01))
summary(modINLA)



MC_samples <- inla.iidkd.sample(10^4, modINLA, "IDrandominterceptLongi", return.cov=TRUE) # for random-effects covariance terms
VarCov <- matrix(unlist(MC_samples), nrow = 2^2)
VarCovMeans <- matrix(rowMeans(VarCov),2,2)
VarCovSD <- matrix(apply(VarCov,1,sd),2,2)
VarCov025 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.025)),2,2)
VarCov975 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.975)),2,2)


















