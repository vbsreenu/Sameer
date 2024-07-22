#!/usr/bin/env Rscript
args <- (commandArgs(TRUE))



ent <- read.table(args[1])
len <- max(ent$V1)

jpeg(file="entropy.jpg", height=400, width=1300, pointsize = 25)
par(mar=c(4,4,1,0))

plot(ent$V1,ent$V2, type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,2),bty="n",lwd=1, main="Entropy of the reads")
axis(2, lwd = 2, col="red", cex.axis=.75)
axis(1,at=seq(1, len, by=(floor(len/11))), col="gray",lwd = 2, cex.axis=.75)


mtext("Genome Length",side=1,line=2.5)
mtext("Entropy",side=2,line=2.5)

# Fush output to PDF
dev.off()
