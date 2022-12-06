library(mvtnorm)
sigma = matrix(rnorm(16,sd = 100),4,4)
sigma = sigma%*%t(sigma)
sigma = sigma/outer(sqrt(diag(sigma)),sqrt(diag(sigma)))

y = rmvnorm(n = 1e7,mean = rep(0,4),sigma = sigma)
mean(y[,1]+sin(y[,2])+exp(y[,3])+y[,4]^2)

x = fastGHQuad::gaussHermiteData(9)$x
w = fastGHQuad::gaussHermiteData(9)$w

X = as.matrix(expand.grid(x,x,x,x))
W = expand.grid(w,w,w,w)
WW = W[,1]*W[,2]*W[,3]*W[,4]
L = chol(sigma)
Lx = sqrt(2)*X%*%L
1/pi^2*sum(WW*(Lx[,1]+sin(Lx[,2])+ exp(Lx[,3])+Lx[,4]^2))

0 + sin(0) + exp(0+0)
