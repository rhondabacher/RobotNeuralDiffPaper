setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

library(Trendy)

# Dynamic peak genes:
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

peak_genes.m <- extractPattern(seg.mouse, Pattern=c("up", "down"), adjR2Cut =.5)
peak_genes.h <- extractPattern(seg.human, Pattern=c("up", "down"), adjR2Cut =.5)

peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

X <- (na.omit(c(peak_genes.h$BreakPoint1)))
Y <- (na.omit(c(peak_genes.m$BreakPoint1)))

wilcox.test(X,Y)
PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
  PP
library("yarrr")

X <- data.frame( Breakpoint = X, Species = "Human")
Y <- data.frame( Breakpoint = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_PeakTime_Figure3.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, 1, 0))
pirateplot(formula = Breakpoint ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "Minute", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()

X <- (na.omit(c(peak_genes.h$BreakPoint1)))
Y <- (na.omit(c(peak_genes.m$BreakPoint1)))

pdf("PLOTS/histogram_peakTimes_Both.pdf", height=6, width=9)
par(mar=c(3,5,1,1), mfrow=c(2,1))
hist(Y, xlim=c(0,600), ylim=c(0,40), 
	col="cornflowerblue", breaks = seq(0, 600, length.out=100), 
	ylab="", main="", xlab="",border="dodgerblue3", 
	cex.axis=1.5, cex.lab=2)
hist(X, xlim=c(0,600), ylim=c(0,15), 
	col="brown1", breaks = seq(0, 600, length.out=100), 
	ylab="", main="", xlab="",border="brown3", 
	cex.axis=1.5, cex.lab=2)
dev.off()


## Spectrum Plots
X <- round(peak_genes.h$BreakPoint1)
Y <- round(peak_genes.m$BreakPoint1)

pdf("PLOTS/spectrum_peakTimes_Both.pdf", height=3, width=13)
par(mar=c(3,3,1,1), mfrow=c(1,1))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', ylab="", xlab="",main="", cex.axis=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=2, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=2, col="cornflowerblue")
dev.off()


## Tables for supplement:

tosave <- formatResults(res.top.h)
tosave <- tosave[peak_genes.h$Gene, ] #order by time of first breakpoint
write.table(tosave, file="TABLES/HumanPeaks_Fig3Genes.csv", row.names=F, quote=F, sep=",")


tosave <- formatResults(res.top.m)
tosave <- tosave[peak_genes.m$Gene, ] #order by time of first breakpoint
write.table(tosave, file="TABLES/MousePeaks_Fig3Genes.csv", row.names=F, quote=F, sep=",")
