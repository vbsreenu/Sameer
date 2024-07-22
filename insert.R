#!/usr/bin/env Rscript
args <- (commandArgs(TRUE))



cov <- read.table(args[1])

c_min <- min(cov$V2)
c_max <- max(cov$V2)
c_max <- c_max+floor((c_max/10))

jpeg(file="insert.jpg", height=400, width=1300, pointsize = 25)
par(mar=c(4,4,1,1))

plot(cov$V1,cov$V2,type="l",col="blue",xlab="",ylab="",ylim=c(c_min,c_max),bty="n",lwd=2, xaxs="i",yaxs="i", xlim=c(c_min,1000))

mtext("Number of pairs",side=2,line=2.5)
mtext("Insert Length",side=1,line=2.5)

# Fush output to PDF
dev.off()
