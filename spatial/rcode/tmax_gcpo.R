### load script to build the mesh
system.time(source('us_mesh.R'))

### load script to get the data
system.time(source('us_get_tmax1day.R'))

### summary of the data
ls()
summary(tmax1day)
sd(tmax1day$tmax, na.rm=TRUE)

### build a mesh

### Construct latent model components
matern <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, 
    prior.sigma = c(5, 0.01), 
    prior.range = c(25, 0.01))

### prior for the likelihood parameter 
lik.prec <- list(prec=list(prior='pc.prec', param=c(5, 0.01)))

### set some INLA parameters
inla.setOption(
    inla.mode='experimental',
    num.threads='4:-1',
    smtp='pardiso',
    inla.call='remote',
    pardiso.license='~/.pardiso.lic')

### projector matrix
As <- inla.spde.make.A(
    mesh, coordinates(tmax1day))

### data stack
stack <- inla.stack(
    tag='data',
    data=list(y=tmax1day$tmax),
    effects=list(list(b0=1, s=1:matern$n.spde)),
    A=list(As))

### fit the model
fit <- inla(
    y ~ 0+b0+f(s, model=matern),
    data=inla.stack.data(stack),
    control.family=list(hyper=lik.prec),
    control.predictor=list(
        A=inla.stack.A(stack),
        compute=TRUE),
    verbose=!TRUE,
    inla.call='remote',
    control.compute=list(cpo=TRUE, po=TRUE))

fit$cpu

### spatial posterior mean
s.mean <- fit$summary.ran$s$mean

### project it into a grid (for plotting) 
y.m <- inla.mesh.project(grid.proj, field=s.mean)
y.m[id.grid.out] <- NA

library(fields)

### visualize the random field + b0
par(mfrow=c(1,1), mar=c(0,0,0,0))
image.plot(
    x=grid.proj$x,
    y=grid.proj$y,
    z=y.m+fit$summary.fix$mean[1], asp=1)
points(tmax1day, cex=0.05, pch=8)
plot(map.moll, add=TRUE, border=gray(0.3,0.5))

### gcpo
system.time(gcpo5 <- inla.group.cv(fit, 5))
system.time(gcpo20 <- inla.group.cv(fit, 20))

### the sum
c(-sum(log(fit$po$po), na.rm=TRUE),
  -sum(log(fit$cpo$cpo), na.rm=TRUE),
  -sum(log(gcpo5$cv), na.rm=TRUE),
  -sum(log(gcpo20$cv), na.rm=TRUE))

summary(sapply(gcpo5$groups, function(x) length(x$idx)))
summary(sapply(gcpo20$groups, function(x) length(x$idx)))

### plot the neighbors for some data points
isel <- c(774, 992, 1714, 2376, 2520, 3252, 3338, 4378, 4596, 5103)

### number of neighbors (with m=10) at the selected data locations
sapply(gcpo20$groups[isel], function(x) length(x$idx)-1)
ce <- coordinates(tmax1day)

par(mfrow=c(1,1), mar=c(0,0,0,0))
image.plot(
    x=grid.proj$x,
    y=grid.proj$y,
    z=y.m+fit$summary.fix$mean[1], asp=1)
plot(map.moll, add=TRUE, border=gray(0.3,0.5))
points(tmax1day, cex=0.5, pch=8)
for(i in isel)
    segments(ce[i, 1], ce[i, 2],
             ce[gcpo20$groups[[i]]$idx[-1], 1],
             ce[gcpo20$groups[[i]]$idx[-1], 2])
points(tmax1day[isel, ], pch=19, cex=2, col='white')


