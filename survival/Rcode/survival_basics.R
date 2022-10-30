scale = 1
shape = 1.5
time_start = 0
time = seq(time_start,10,length.out = 1000)
pdf = dnorm(time)/pnorm(0,lower.tail = F)
survival = pnorm(time,lower.tail = FALSE)/pnorm(0,lower.tail = F)
# pdf = dexp(time)
# survival = pexp(time,lower.tail = FALSE)

harzard = pdf/survival
cumulative_harzard = unlist(lapply(1:length(time),FUN = function(i){pracma::trapz(time[1:i],harzard[1:i])}))
plot(time,pdf,type="l")
plot(time,survival,type="l")
plot(time,harzard,type="l")
plot(time,cumulative_harzard,type = "l")

start_time = 4
time = seq(time_start,10,length.out = 1000)
