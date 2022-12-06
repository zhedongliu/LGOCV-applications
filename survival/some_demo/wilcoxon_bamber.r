n = 10000
eta = rnorm(n,sd = 1)
p = 1/(1+exp(-eta))
y = rbinom(n = n,size = 1,prob = p)
TPR = function(c){
    sum(y*(p>c))/sum(y)
}
FPR = function(c){
    sum((1-y)*(p>c))/sum((1-y))
}
c = seq(1,0,-0.0001)
xx = unlist(lapply(c,FPR))
yy = unlist(lapply(c,TPR))
plot(xx,yy,type="l")
abline(a=0,b=1)
# AUC = pracma::trapz(unlist(lapply(c,FPR)),unlist(lapply(c,TPR)))

true_fold = which(y == 1)
false_fold = which(y == 0)
W = 0
for(i in true_fold){
    W = W + sum((p[false_fold]<p[i]))
}
AUC_wilcoxon = W/sum(y)/sum(1-y)
AUC_wilcoxon
# AUC


# upper = 0
# lower = 0
# for(i in 1:n){
#     for(j in 1:n){
#         upper = upper + (p[i] > p[j])*(y[i] == 1)*(1- y[j] == 1)
#         lower = lower + (y[i] == 1)*(1- y[j] == 1)
#     }
# }
# upper/lower
# AUC_wilcoxon
# AUC
