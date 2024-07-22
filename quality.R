#!/usr/bin/env Rscript
args <- (commandArgs(TRUE))

qual <- read.table(args[1])

mapped_max <- max(qual$V2)
mapped_min <- min(qual$V2)

unmapped_max <- max(qual$V3)
unmapped_min <- min(qual$V3)

jpeg(file="quality.jpg", height=400, width=1300, pointsize = 25)
par(mar=c(4,4,1,4))

plot(qual$V1,qual$V2,type="l",col="blue",yaxt="n",xlab="",ylab="",ylim=c(mapped_min,mapped_max),bty="n",lwd=2, main="Read quality")
axis(2, col="blue",lwd = 3)
mtext("Number of mapped Reads",side=2,line=2.5)

par(new=TRUE)

plot(qual$V3,type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(unmapped_min,unmapped_max),bty="n",lwd=1)
axis(4, col="red",lwd = 3)
mtext("Number of unmapped Reads",side=4,line=2.5)

#axis(1,qual$V1,col="gray",lwd = 2)
mtext("Read Quality",side=1,line=2.5)

# Fush output to PDF
dev.off()
