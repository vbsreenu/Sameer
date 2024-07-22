#!/usr/bin/env Rscript
args <- (commandArgs(TRUE))

map <- read.table(args[1])
map_max <- max(map$V2)

map_max=map_max + (map_max/10)

mapped=head(map$V2,1); unmapped=tail(head(map$V2,2),1)
map_per=round((mapped/(mapped+unmapped))*100)
unmap_per=round((unmapped/(mapped+unmapped))*100)

jpeg(file="map.jpg", height=300, width=250, pointsize = 12)
par(mar=c(2,4,1,0))

barplot(c(map_per,unmap_per),  col=c("blue","red"), names.arg=c("Mapped", "Unmapped"), ylim=c(0,100), space=.75,ylab="Percentage", main="Read Composition")

# Fush output to PDF
dev.off()
