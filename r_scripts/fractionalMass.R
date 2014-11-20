#options
options(digits=10)

if(length(commandArgs(TRUE)) == 0){
    write(R.version$version.string, stderr())
    write("script version: no version available. SVN revision number will be included.",stderr())
} else{
    csvin<-commandArgs(TRUE)[1]
    out<-commandArgs(TRUE)[2]
    theo_mass <-commandArgs(TRUE)[3]
    #post<-commandArgs(TRUE)[2]

    R <- read.delim(csvin, header = TRUE);
    png(out,width=480, height=320, pointsize=12, bg="#FFFFFF", res=NA);
    theoretical_masses_var<-paste("file",theo_mass,sep=":")
    theoretical_masses_file<-gsub("file:","", theoretical_masses_var)
    b<-read.delim(theoretical_masses_file)
    a<-R
    a<-a$MZ*a$Charge - a$Charge*1.007276467
    a1<-floor(a);
    a2<-a-floor(a);
    a<-cbind(a1,a2);
    smoothScatter(b,nrpoints=0,ylab="fractional mass",xlab="nominal mass", xlim=c(0,5000),pch=4,col="blue")
    points(a,col="#88000011", pch=20,cex=1)
    legend("topleft",c("theoretical","experimental"),col=c("darkblue","darkred"),pch=19,bty='n')
    dev.off();
}
