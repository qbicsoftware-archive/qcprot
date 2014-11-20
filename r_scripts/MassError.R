#options
options(digits=10)

csvin<-commandArgs(TRUE)[1]
out<-commandArgs(TRUE)[2]
#post<-commandArgs(TRUE)[2]



R <- read.delim(csvin, header = TRUE);
png(out,width=480, height=320, pointsize=12, bg="#FFFFFF", res=NA);

plot(R$"RT"/60,R$"delta_ppm",xlab="RT (min)",ylab="mass error in [ppm]",main="", ylim=c(-11,11), pch='.')
lines(lowess(R$"RT"/60,R$"delta_ppm", f=1/20), col="red")
abline(h=0, col="blue")

dev.off();