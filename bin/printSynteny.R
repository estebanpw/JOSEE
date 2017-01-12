##!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# library(lattice)
# 
# if(length(args) < 3){
#   stop("USE: Rscript --vanilla fragsTrimPlot.R <results_absolute> <n_files> <distance>")
# }
# 
# 
# path <- args[1]
# anexo <- ".trim.csv_"
# n_files <- as.numeric(args[2])
# current <- 1
# open <- 0
# close <- 0
# distance <- as.numeric(args[3])
# first_time <- 1
# line_width = 10
# max_len <- 1000000
# pathOut <- paste(args[1], paste("_d_", paste(distance, ".png", sep=""), sep=""), sep="")
# 
# png(pathOut, width = 800, height = 400)

path <- "C:/Users/Annette/Documents/Cosas de esteban/Jose/data.csv"
synteny <- scan(path, sep="\n", what=character())


# 1 -> start
# 2 -> end
# 3 -> order
# 4 -> length
# 5 -> genome

curr_blocks_per_synteny <- 0
blocks <- c()
current_block <- 1
mycols = c("red", "black", "green", "blue", "yellow")
total <- 1
was_se_next = 0

for(k in 1:length(synteny)){
  
  if(synteny[k] == "NEXT" && was_set_next == 0){
    was_set_next = 1
  }
  if(synteny[k] == "BREAK"){
    was_set_next = 0
    
    # Find min and max
    c_min <- 10000000
    c_max <- 0
    for(j in 1:(current_block-1)){

      if(c_min > as.numeric(blocks[[j]][1])){
        c_min <- as.numeric(blocks[[j]][1])
      }
    }
    for(j in 1:(current_block-1)){
      if(c_max < as.numeric(blocks[[j]][2])){
        c_max <- as.numeric(blocks[[j]][2])
      }
    }
    
    
    # Draw time
    png(paste(paste("out_", total, sep=""), ".png", sep=""), width = 800, height = 400)
    plot(c(c_min, c_max), c(1, 100), type= "n", xlab = "", ylab = "")
    # rect(xleft, ybottom, xright, ytop
    curr_g <- 0
    for(j in 1:(current_block-1)){
      start <- as.numeric(blocks[[j]][1])
      end <- as.numeric(blocks[[j]][2])
      pos_low <- as.numeric(blocks[[j]][5])*10
      pos_high <- (as.numeric(blocks[[j]][5])+1)*10 
      
      if((j-1) %% n_genomes == 0){
        curr_g <- curr_g + 1
      }
      rect(start, pos_low, end, pos_high, col = mycols[curr_g])
    }
    total <- total + 1
    dev.off()
    # Reset blocks
    blocks <- c()
    current_block <- 1
  }else{
    
    blocks[current_block] <- strsplit(synteny[k], split = ";")
    current_block <- current_block + 1
  }
  
  
}

# for(k in 0:n_files){
#   current <- 1
#   
#   profile <- scan(paste(path, paste(anexo, k, sep=""), sep=""), sep="\n", what=character())
#   
#   if(first_time == 1){
#     plot(c(1, length(profile)*60), c(0,0), col = "black", ylim=c(0, n_files/2),  xlim=c(1, max_len), lwd=line_width, type="l", xlab="Genome", ylab="Blocks")
#     first_time = 0
#   }else{
#     lines(c(1, length(profile)*60), c(k/2,k/2), lwd=line_width,  col="black")  
#   }
#   
#   
#   for(i in 1:length(profile)){
#     cs = strsplit(profile[i], split = character(0))
#     
#     #c holds each character
#     
#     
#     for(j in 1:length(cs[[1]])){
#       current = current + 1
#       
#       if(cs[[1]][j] == "1"){
#         open = current
#       } 
#       if(cs[[1]][j] == "3"){
#         close = current
#       } 
#       
#       if((close - open) > distance){
#         lines(c(open, current), c(k/2,k/2), lwd=2,  col="red")
#         
#       }
#       
#     }
#     
#   }
# }
# dev.off()
# 


# 
# png(pathOut, width = xsize, height = ysize)
# plot(1:points*swsize, array, type="l", xlab= xlabel, ylab= ylabel, main="Genome profiling of the accumulated number of nucleotides per base")
# dev.off()

