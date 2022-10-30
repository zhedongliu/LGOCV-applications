library(spam)
t = seq(0,10,0.01)
n = length(t)
#hazard = abs(rnorm(length(t),mean = t))
hazard = exp(sin(t))

cumulative_harzard = numeric(length(t))
for(i in 2:length(t)){
    cumulative_harzard[i] = cumulative_harzard[i-1]+hazard[i]*0.01
}
survival = exp(-cumulative_harzard)
pdf = hazard*survival
plot(t,pdf,type="l")
plot(t,survival,type="l")
pracma::trapz(t,pdf)
