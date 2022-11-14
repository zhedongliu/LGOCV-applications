n = 10000
eta = rnorm(n)
p = 1/(1+exp(-eta))
y = rbinom(n = n,size = 1,prob = p)
TPR = function(c){
    sum(y*(p>c))/sum(y)
}
FPR = function(c){
    sum((1-y)*(p>c))/sum((1-y))
}

c = seq(1,0,-0.001)
plot(unlist(lapply(c,FPR)),unlist(lapply(c,TPR)),type="l")
AUC = pracma::trapz(unlist(lapply(c,FPR)),unlist(lapply(c,TPR)))

true_fold = which(y == 1)
false_fold = which(y == 0)
W = 0
for(i in true_fold){
    W = W + sum((p[false_fold]<p[i]))
}
AUC_wilcoxon = W/sum(y)/sum(1-y)
AUC
AUC_wilcoxon