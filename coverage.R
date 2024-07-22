#!/usr/bin/env Rscript
args <- (commandArgs(TRUE))



cov <- read.table(args[1])

c_max <- max(c(max(cov$V2),abs(min(cov$V3))))
c_max <- c_max+floor((c_max/10))
len <- max(cov$V1)

jpeg(file="coverage.jpg", height=400, width=1300, pointsize = 25)
par(mar=c(4,4,1,0))

plot(cov$V1,cov$V2,type="l",xaxt="n",yaxt="n",col="blue",xlab="",ylab="",ylim=c(-c_max,c_max),bty="n",lwd=2, main="Coverage by reads")

yt<-seq(-c_max,0,by=floor((c_max/2)))
axis(2, col="gray",lwd = 2,cex.axis=.5,at=yt)

yt<-seq(0,c_max,by=floor((c_max/2)))
axis(2, col="blue",lwd = 2,cex.axis=.5,at=yt)

mtext(c("+ve reads","-ve reads"),side=2,line=2.5,at=c((c_max/1.75),-(c_max/1.75)),col=c("blue","azure4"))


abline(h = 0, col = "red", lty=3)

par(new=TRUE)


plot(cov$V1,cov$V3,type="l",col="gray",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-c_max,c_max),bty="n",lwd=2)

axis(1,at=seq(1, len, by=(floor(len/11))), col="gray",lwd = 2, cex.axis=.75)

mtext("Genome Length",side=1,line=2.5)

# Fush output to PDF
dev.off()
