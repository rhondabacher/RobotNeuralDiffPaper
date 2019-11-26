setwd("~/RobotNeuralDiffPaper/")


load("RDATA/jointPlots_loadDataBoth.Rdata")




## First get slopes Up and Down for ANY PEAK GENE:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

top.mouse <- top.mouse[peak_genes.m[,1],]
top.human <- top.human[peak_genes.h[,1],]

peak.com.h <- peak_genes.h
whichSlope.up.h<-c()
whichSlope.down.h<-c()
for(i in 1:nrow(peak.com.h)) {
	keep <- which(round(res.top.h$Breakpoints[peak.com.h[i,1],]) == round(c(peak.com.h[i,2])))
	whichSlope.up.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],keep]
	whichSlope.down.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],(keep+1)]
	
}
names(whichSlope.up.h) <- peak.com.h[,1]
names(whichSlope.down.h) <- peak.com.h[,1]


for(i in names(which(whichSlope.down.h > 0))) { 
	keep <- which(round(res.top.h$Breakpoints[i,]) == round(c(peak.com.h[which(peak.com.h[,1] == i),2])))
	whichSlope.down.h[i] <- res.top.h$Segment.Slopes[i,(keep+2)]
}

peak.com.m <- peak_genes.m
whichSlope.up.m<-c()
whichSlope.down.m<-c()
for(i in 1:nrow(peak.com.m)) {
	keep <- which(round(res.top.m$Breakpoints[peak.com.m[i,1],]) == round(c(peak.com.m[i,2])))
	whichSlope.up.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],keep]
	whichSlope.down.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],(keep+1)]	
}
names(whichSlope.up.m) <- peak.com.m[,1]
names(whichSlope.down.m) <- peak.com.m[,1]

for(i in names(which(whichSlope.down.m > 0))) { 
	keep <- which(round(res.top.m$Breakpoints[i,]) == round(c(peak.com.m[which(peak.com.m[,1] == i),2])))
	whichSlope.down.m[i] <- res.top.m$Segment.Slopes[i,(keep+2)]
}

X <- whichSlope.up.h
Y <- whichSlope.up.m
mean(X)
mean(Y)

# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = abs(whichSlope.up.h), Species = "Human")
X <- data.frame( Slope = abs(whichSlope.up.m), Species = "Mouse")

longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals<- pairwiseCI(Slope ~ Species,
           data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

# pdf("PLOTS/histogram_allPeak_slopeUP_Figure6.pdf", height=6, width=7.5)
# par(mar=c(6,6,3,1))
# hist(X, xlim=c(0,.05), ylim=c(0,600), border="brown3",
# 	col=alpha("brown1", .6), breaks = seq(0, .05, length.out=50),
# 	main="", xlab="Slope (Scaled Expression/Minute)",
# 	cex.axis=1.8, cex.lab=2)
# hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, .053, length.out=50))
# legend('topright', c("Mouse","Human"), lwd=3, col=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2)
# dev.off()


library("yarrr")

# X <- data.frame( Slope = whichSlope.up.h, Species = "Human")
# Y <- data.frame( Slope = whichSlope.up.m, Species = "Mouse")
#
# longdata <- rbind(Y, X)
X <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_allPeak_slopeUP_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()



X <- whichSlope.down.h
Y <- whichSlope.down.m
mean(X)
mean(Y)


# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
X <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")

longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)


X <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_allPeak_slopeDown_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()




########################################################################################################
########################################################################################################
########################################################################################################

## Now do this for Ortholog peaks:
								
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)

# Get Up and Down Slope:
peak.com.h <- subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)
whichSlope.up.h<-c()
whichSlope.down.h<-c()
for(i in 1:nrow(peak.com.h)) {
	keep <- which(round(res.top.h$Breakpoints[peak.com.h[i,1],]) == round(c(peak.com.h[i,2])))
	whichSlope.up.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],keep]
	whichSlope.down.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],(keep+1)]
	
}
names(whichSlope.up.h) <- peak.com.h[,1]
names(whichSlope.down.h) <- peak.com.h[,1]

peak.com.m <- subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)
whichSlope.up.m<-c()
whichSlope.down.m<-c()
for(i in 1:nrow(peak.com.m)) {
	keep <- which(round(res.top.m$Breakpoints[peak.com.m[i,1],]) == round(c(peak.com.m[i,2])))
	whichSlope.up.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],keep]
	whichSlope.down.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],(keep+1)]	
}
names(whichSlope.up.m) <- peak.com.m[,1]
names(whichSlope.down.m) <- peak.com.m[,1]



for(i in names(which(whichSlope.down.h > 0))) { 
	keep <- which(round(res.top.h$Breakpoints[i,]) == round(c(peak.com.h[which(peak.com.h[,1] == i),2])))
	whichSlope.down.h[i] <- res.top.h$Segment.Slopes[i,(keep+2)]
}
for(i in names(which(whichSlope.down.m > 0))) { 
	keep <- which(round(res.top.m$Breakpoints[i,]) == round(c(peak.com.m[which(peak.com.m[,1] == i),2])))
	whichSlope.down.m[i] <- res.top.m$Segment.Slopes[i,(keep+2)]
}


X <- whichSlope.up.h
Y <- whichSlope.up.m

# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
X <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")

longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

X <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_commonPeak_slopeUp_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()





X <- whichSlope.down.h
Y <- whichSlope.down.m

# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
X <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")

longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

X <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_commonPeak_slopeDown_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()




#######################################################################################################
########################################################################################################
########################################################################################################

## Now do this for all Orthologs:

library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

all.slopes.mouse <- res.top.m$Segment.Slopes
colnames(all.slopes.mouse) <- paste0("mouse_slope", seq_len(ncol(all.slopes.mouse)))

all.slopes.human <- res.top.h$Segment.Slopes
colnames(all.slopes.human) <- paste0("human_slope", seq_len(ncol(all.slopes.human)))


all.slopes.human <- (all.slopes.human[ortho.genes.use$hgnc_symbol,1])
all.slopes.mouse <- (all.slopes.mouse[ortho.genes.use$mgi_symbol,1])

all.slopes.human <- all.slopes.human[names(which(res.top.h$Segment.Trends[names(all.slopes.human),1] !=0))]
all.slopes.mouse <- all.slopes.mouse[names(which(res.top.m$Segment.Trends[names(all.slopes.mouse),1] !=0))]

whichSlope.up.h <- all.slopes.human[all.slopes.human>0]
whichSlope.up.m <- all.slopes.mouse[all.slopes.mouse>0]

# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
X <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")

longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

X <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_commonGenes_slopeUp_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()






whichSlope.down.h <- all.slopes.human[all.slopes.human<0]
whichSlope.down.m <- all.slopes.mouse[all.slopes.mouse<0]

# Slope meaning is not easy to understand. Use relative effect:
Y <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
X <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")
longdata <- rbind(Y, X)

library(pairwiseCI)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

X <- data.frame( Slope = (whichSlope.down.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.down.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_commonGenes_slopeDown_Figure6.pdf", height=2, width=2, useDingbats=F, family = "Arial")
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()


