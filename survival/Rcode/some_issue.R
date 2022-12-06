y <- rep(NA, 15)
y[1:5] <- rpois(5, lambda = 1)

yy <- rep(NA, 15)
yy[6:10] <- rnorm(5)

yyy <- inla.surv(c(rep(NA, 10), 1:5), c(rep(NA, 10), rep(1, 5)), cure = cbind(1, 1:15))

r <- inla(Y ~ 1,
          family = c("exponential", "normal", "weibullsurv"), 
          data = list(Y = list(y = y, yy = yy, yyy = yyy)), 
          control.fixed = list(prec.intercept = 1), 
          control.family = list(list(), list(), list(variant = 1,hyper = list(theta1 = list(initial = 2,fixed = TRUE)))), 
          ## will enable config=T
          control.compute = list(likelihood.info = TRUE),
          verbose = FALSE,
          safe = FALSE)

#linear.predictor should be outside of a11$family.arg.str 
a11=eval(parse(text=r$misc$configs$config[[1]]$arg.str[11]))
a11$linear.predictor = a11$family.arg.str$linear.predictor
a11$y.surv$lower = a11$y.surv$time
a11$y.surv$time = 0
a11$y.surv$event = 0
a11$cure.prob = 0
do.call(inla.likelihood, args = a11)
pweibull(q = a11$y.surv$lower,shape = a11$family.arg.str$alpha,scale = exp(a11$family.arg.str$linear.predictor),lower.tail = FALSE,log.p = T)

r$misc$configs$config[[1]]$arg.str[11]
r$misc$configs$config[[1]]
r$summary.hyperpar
