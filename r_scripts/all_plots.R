#!/bin/env Rscript

# Make quality control plots from csv files extracted from a qcML file.
#
# Usage: Rscript all_plots.R <csvdir> <outdir>
#
# <indir> must contain the following csv files:
#     mass_accuracy.csv :
#         Extracted from accession QC:0000038. Must contain the colums
#         ``MZ`` and ``Charge``


options(digits=10)

#-------------------------Fractional Masses
fractional_masses <- function(csv_dir, out_dir) {
  print(dir(csv_dir))
  
  # get theoretical masses from the directory in which this script is stored
  argv <- commandArgs(trailingOnly = FALSE)
  script_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  theo_mass_path <- paste(script_dir, "theoretical_masses.txt", sep="/")
  theo_mass <- read.delim(theo_mass_path)

  mass_accuracy <- read.delim(file.path(csv_dir, "features.csv"),
                              header = TRUE)

  # make plot
  name = "01_fractional_masses"
  mass <- (mass_accuracy$MZ * mass_accuracy$Charge
           - mass_accuracy$Charge * 1.007276467)
  integer_part <- floor(mass)
  fraction <- mass - integer_part

  df <- cbind(integer_part, fraction)

  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width = 640, height = 640, pointsize = 12, bg = "#FFFFFF", res = NA)

  smoothScatter(theo_mass, nrpoints = 0, ylab = "fractional mass",
                 xlab = "nominal mass", xlim = c(0,5000), pch = 4, col = "darkblue")
  points(df, col="darkred", pch = 20, cex = 1)
  #plot(df[,1], df[,2],   pch = 20, cex = 1)
  legend("topleft", c("theoretical", "experimental"),
         col = c("darkblue", "darkred"), pch = 19, bty = 'n')
  dev.off()

  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Fractional masses")
  sink()

  # write description
}

##---------------------MASS ACCURACY
mass_accuracy <- function(csv_dir, out_dir) {
  mass_accuracy <- read.delim(file.path(csv_dir, "identifications.csv"),header = TRUE)
  # make plot
  name = "04_mass_accuracy"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width = 640, height = 640, pointsize = 12, bg = "#FFFFFF", res = NA)
  
    hist(mass_accuracy$delta_ppm,xlim=c(-10,10),breaks=seq(min(mass_accuracy$delta_ppm)-0.01, max(mass_accuracy$delta_ppm)+0.01, 0.01),xlab="ppm",
         main=paste("delta ppm"))
    abline(v=median(mass_accuracy$delta_ppm),col="red", lwd=2)
    mtext(paste("median(accuracy)=",round(median(mass_accuracy$delta_ppm),3)," ppm",sep=""))
  dev.off()

  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Mass accuracy")
  sink()

}# mass_accuracy

#------------------MASS ERROR
mass_error <- function(csv_dir, out_dir) {
  mass_error <- read.delim(file.path(csv_dir, "mass_error.csv"),header = TRUE)

  name = "03_mass_error"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width=640, height=320, pointsize=12, bg="#FFFFFF", res=NA);

  plot(mass_error$"RT"/60,mass_error$"delta_ppm",xlab="RT (min)",ylab="mass error in [ppm]",main="", ylim=c(-11,11), pch='.')
  abline(h=0, col="blue")
  dev.off();
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Mass Error")
  sink()

}# mass_error

#-----------------------TIC

tic <- function(csv_dir, out_dir) {
  tic <- read.delim(file.path(csv_dir, "tic.csv"),header = TRUE)
  colnames(tic) <- cbind("RT","TIC")
  # make plot
  name = "02_tic"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width = 640, height = 640, pointsize = 12, bg = "#FFFFFF", res = NA)
  
  res = barplot(t(tic$TIC), xlab="RT (min)",ylab="Intensity")
  time_seq = seq(min(tic$RT),max(tic$RT),1200)
  time_seq = round(time_seq)/60
  t<-which(time_seq %% 1==0 & duplicated(time_seq)==F)
  tmp = seq(min(res),max(res),max(res)/(length(time_seq)))
  axis(1,at=tmp[t],labels=time_seq[t])
  
  dev.off()
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Total ion current chromatogram (TIC)")
  sink()
  
}# TIC

#------------------Charge Distribution
charge <- function(csv_dir, out_dir) {
  charge <- read.delim(file.path(csv_dir, "features.csv"),header = TRUE)
  name = "05_charge"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width=320, height=320, pointsize=12, bg="#FFFFFF", res=NA);
    
  hist(charge$"Charge", xlab="charge", main="", col="cyan4")
   
  dev.off();
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Charge Distribution of Identified Features")
  sink()
  
}# charge

#------------------Id ratios = identified features vs observed
id_ratios <- function(csv_dir, out_dir) {
  id_ratios <- read.delim(file.path(csv_dir, "features.csv"),header = TRUE)
  #print(head(id_ratios))
  
  id_ratios_raw <- read.delim(file.path(csv_dir, "raw_MS1.csv"),header = TRUE)
  colnames(id_ratios_raw) <- c("RT","Precursor")
  #print(head(id_ratios_raw))
  
  name = "06_id_ratios"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width=640, height=640, pointsize=12, bg="#FFFFFF", res=NA);
  
  plot(id_ratios_raw$RT/60,id_ratios_raw$Precursor,pch=16,,xlab="RT (min)",ylab="m/z",cex=0.3)
  points(id_ratios$RT/60,id_ratios$MZ,col="red",pch=4,cex=0.3)
 
  legend("topleft",c("Recorded spectra","Identified spectra"),pch=19,col=c(1,2))
   
  
  dev.off();
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Identified Features")
  sink()
  
}# id_ratios

#------------------Injection times MS1 and MS2

injection_times <- function(csv_dir, out_dir) {
  
  inj_times <- read.delim(file.path(csv_dir, "injection_times.csv"), header = TRUE, sep=",")

  #ms1
  name = "07_injection_times_ms1"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width=640, height=440, pointsize=12, bg="#FFFFFF", res=NA);
  
  ms1 = subset(inj_times, mslevel==1)
  plot(ms1$rt, ms1$time, xlab="RT (minutes)", ylab="Ion Injection Time (ms)", cex=0.3,  main="MS 1")
  
  dev.off();
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Injection Times MS1")
  sink()
  
  #ms2
  name = "08_injection_times_ms2"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width=640, height=440, pointsize=12, bg="#FFFFFF", res=NA);
  
  ms2 = subset(inj_times, mslevel==2)
  plot(ms2$rt, ms2$time, xlab="RT (minutes)", ylab="Ion Injection Time (ms)", cex=0.3,  main="MS 2")
  
  dev.off();
  # write title
  sink(file.path(out_dir, paste(name, ".title")))
  cat("Injection Times MS2")
  sink()
  
}# end --Injection times
#----------------------------------------


if (length(commandArgs(TRUE)) != 2) {
  write("Invalid arguments. Usage: Rscript all_plots.R <csvdir> <outdir>",
        stderr())
  quit(save = "no", status = 1)
} else {
  
  csv_dir <- commandArgs(TRUE)[1]
  out_dir <- commandArgs(TRUE)[2]
  
  mass_accuracy(csv_dir, out_dir)
  mass_error(csv_dir, out_dir)
  fractional_masses(csv_dir, out_dir)
  tic(csv_dir, out_dir)
  charge(csv_dir, out_dir)
  id_ratios(csv_dir, out_dir)
  injection_times(csv_dir, out_dir)
 
}
