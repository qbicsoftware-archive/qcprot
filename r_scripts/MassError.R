#options
options(digits=10)

csvin<-commandArgs(TRUE)[1]
out<-commandArgs(TRUE)[2]
#post<-commandArgs(TRUE)[2]



R <- read.delim(csvin, header = TRUE);
png(out,width=640, height=320, pointsize=12, bg="#FFFFFF", res=NA);

plot(R$"RT"/60,R$"delta_ppm",xlab="RT (min)",ylab="mass error in [ppm]",main="", ylim=c(-11,11), pch='.')
abline(h=0, col="blue")

dev.off();
