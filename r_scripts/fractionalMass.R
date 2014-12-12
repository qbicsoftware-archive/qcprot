#!/bin/env Rscript

# Make quality control plots from csv files extracted from a qcML file.
#
# Usage: Rscript fractionalMass.R <indir> <outdir>
#
# <indir> must contain the following csv files:
#     mass_accuracy.csv :
#         Extracted from accession QC:0000038. Must contain the colums
#         ``MZ`` and ``Charge``

options(digits=10)

fractional_masses_plot <- function(csv_dir, out_dir) {
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
  png(plot_path, width = 640, height = 640, pointsize = 12,
      bg = "#FFFFFF", res = NA)

  smoothScatter(theo_mass, nrpoints = 0, ylab = "fractional mass",
                 xlab = "nominal mass", xlim = c(0,5000), pch = 4, col = "blue")
print(head(df))
print(dim(df))
  points(df, col = "red", pch = 20, cex = 1)
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

##---------------------
mass_accuracy <- function(csv_dir, out_dir) {
  mass_accuracy <- read.delim(file.path(csv_dir, "identifications.csv"),
                              header = TRUE)
  print(head(mass_accuracy))
  # make plot
  name = "02_mass_accuracy"
  plot_path = file.path(out_dir, paste(name, ".png"))
  png(plot_path, width = 640, height = 640, pointsize = 12,
      bg = "#FFFFFF", res = NA)
  
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

if (length(commandArgs(TRUE)) != 2) {
  write("Invalid arguments. Usage: Rscript fractionalMass.R <indir> <outdir>",
        stderr())
  quit(save = "no", status = 1)
} else {
  csv_dir <- commandArgs(TRUE)[1]
  out_dir <- commandArgs(TRUE)[2]

  fractional_masses_plot(csv_dir, out_dir)
  mass_accuracy(csv_dir, out_dir)
}
