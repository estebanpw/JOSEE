#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(lattice)

if(length(args) < 3){
  stop("USE: Rscript --vanilla fragsTrimPlot.R <results_absolute> <n_files> <distance>")
}


path <- args[1]
anexo <- ".trim.csv_"
n_files <- as.numeric(args[2])
current <- 1
open <- 0
close <- 0
distance <- as.numeric(args[3])
first_time <- 1
line_width = 10
max_len <- 1000000
pathOut <- paste(args[1], paste("_d_", paste(distance, ".png", sep=""), sep=""), sep="")

png(pathOut, width = 800, height = 400)

for(k in 0:n_files){
  current <- 1
  
  profile <- scan(paste(path, paste(anexo, k, sep=""), sep=""), sep="\n", what=character())
  
  if(first_time == 1){
    plot(c(1, length(profile)*60), c(0,0), col = "black", ylim=c(0, n_files/2),  xlim=c(1, max_len), lwd=line_width, type="l", xlab="Genome", ylab="Blocks")
    first_time = 0
  }else{
    lines(c(1, length(profile)*60), c(k/2,k/2), lwd=line_width,  col="black")  
  }
  
  
  for(i in 1:length(profile)){
    cs = strsplit(profile[i], split = character(0))
    
    #c holds each character
    
    
    for(j in 1:length(cs[[1]])){
      current = current + 1
      
      if(cs[[1]][j] == "1"){
        open = current
      } 
      if(cs[[1]][j] == "3"){
        close = current
      } 
      
      if((close - open) > distance){
        lines(c(open, current), c(k/2,k/2), lwd=2,  col="red")
        
      }
      
    }
    
  }
}
dev.off()



# 
# png(pathOut, width = xsize, height = ysize)
# plot(1:points*swsize, array, type="l", xlab= xlabel, ylab= ylabel, main="Genome profiling of the accumulated number of nucleotides per base")
# dev.off()
