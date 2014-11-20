#options
options(digits=10)

csvin_p<-commandArgs(TRUE)[1]

csvin_id<-commandArgs(TRUE)[2]
post<-commandArgs(TRUE)[3]

#MS:1000040 = Precursor m/z
#MS:1000894_[sec] = RT

file_p <- read.delim(csvin_p, header = TRUE);
file_id <- read.delim(csvin_id, header = TRUE);

colnames(file_p) <- cbind("RT","Precursor")

#print(file_p$Precursor)
#R <- read.csv("/tmp/R-inDataTempFile-7432835522670914453.csv", header = TRUE, row.names = 1);
png(post, width=480, height=320, pointsize=12, bg="#FFFFFF", res=NA);
plot(file_p$RT/60,file_p$Precursor,pch=16,,xlab="RT (min)",ylab="m/z",cex=0.3)
points(file_id$RT/60,file_id$MZ,col="red",pch=4,cex=0.3)
legend("topleft",c("recorded spectra","identified spectra"),pch=19,col=c(1,2))

dev.off();