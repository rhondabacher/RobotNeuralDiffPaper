setwd("~/RobotSeq/")

## This code is the analysis of the subset of mouse cells s/t sampling frequency is 12m
## This code generates the figures for Supplementary Figure XX.

load("RDATA/jointPlots_loadDataBoth.Rdata")

rm(list=setdiff(ls(), c("orth.genes.clean", "seg.human", "t.v.h", "data.norm.scale.h", "peak_genes.h")))

library(Trendy)


load("RDATA/trendy_run_mouse_scaleData_shortSubset.RData")
seg.mouse <- results(seg.mouse.scaled.short)

data.norm.scale.m <- data.norm.mouse.short.scale
t.v.m <- t.v.short


library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)

trendy.summary.m <- formatResults(res.top.m)

peak_genes <- extractPattern(seg.mouse, Pattern=c("up", "down"), adjR2Cut =.2)
peak_genes <- rbind(peak_genes, extractPattern(seg.mouse, Pattern=c("up", "same", "down"), adjR2Cut =.2)[,1:2])
peak_genes.m <- peak_genes

## Rerun majority of main plots in analysis for supplementary figure

# For Ortholog Genes, what % of genes have an iniital trend?
# Dynamic genes:
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

top1 <- merge(orth.genes.clean, top.human, by="hgnc_symbol")
top2 <- merge(top1,top.mouse, by="mgi_symbol")
top2 <- top2[!duplicated(top2),][,1:3]

dupg <- top2[which(duplicated(top2[,2])),2] 
subset(top2, hgnc_symbol %in% dupg)
toRM <- c()
for(i in 1:length(dupg)) {

  tempSet <- subset(top2, hgnc_symbol %in% dupg[i])
  bestMatch <- which(tempSet[,2] == toupper(tempSet[,1]))[1]
  if (!is.na(bestMatch)) {
    toRM <- c(toRM, rownames(subset(top2, hgnc_symbol %in% dupg[i]))[-bestMatch])
  }
  else {
    maxOrth <- which.max(subset(top2, hgnc_symbol %in% dupg[i])[,3])
    toRM <- c(toRM, rownames(subset(top2, hgnc_symbol %in% dupg[i]))[-maxOrth])
  }
}
dupg <- top2[which(duplicated(top2[,1])),1] 
subset(top2, mgi_symbol %in% dupg)
for(i in 1:length(dupg)) {

  tempSet <- subset(top2, mgi_symbol %in% dupg[i])
  bestMatch <- which(tempSet[,2] == toupper(tempSet[,1]))[1]
  if (!is.na(bestMatch)) {
    toRM <- c(toRM, rownames(subset(top2, mgi_symbol %in% dupg[i]))[-bestMatch])
  }
  else {
    maxOrth <- which.max(subset(top2, mgi_symbol %in% dupg[i])[,3])
    toRM <- c(toRM, rownames(subset(top2, mgi_symbol %in% dupg[i]))[-maxOrth])
  }
}
toRM <- unique(toRM)

top2 <- top2[-as.numeric(toRM),]
dim(top2) 


# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse, !(Gene %in% top2$mgi_symbol))
mouse.genes.check1 <- mouse.genes.check1[which(toupper(mouse.genes.check1[,1]) %in% top.human$Gene),]
recovered.g <- data.frame(mgi_symbol = mouse.genes.check1[,2], hgnc_symbol = toupper(mouse.genes.check1[,2]), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)

ortho.genes.use <- top2

dim(ortho.genes.use)
# 1112

## All pre-processing above!
########################################################################################################

# Supp Figure A and B and C. On Orthologs!
# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.

timeUpDown.m <- t.v.m[apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) which(x != 0)[1])]

library(ggplot2)

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100

propTestVals <-prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 2)

pdf("PLOTS/percent_FirstTime_UpDown_Orthologs_shortMouse-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(sum(timeUpDown.h == 0), sum(timeUpDown.m==0), sum(timeUpDown.h > 0), sum(timeUpDown.m > 0)),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,1500), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()

## ORTHOLOGS Monotonic genes, only up or down the entire time:

onlyDown.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == 1)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

propTestVals <- prop.test(c((all.mono.h), (all.mono.m)), c(length(ortho.genes.use$mgi_symbol), length(ortho.genes.use$mgi_symbol)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 2)

pdf("PLOTS/numGenes_Monotonic_Orthologs_shortMouse-SuppFig1.pdf", height=2,width=2.5)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,650), cex.axis=.6, cex.lab=.7
)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()

all.mono.h
all.mono.m
all.mono.h / all.mono.m






## Get breakpoints and slopes for ortholog genes

library(gplots)
library(ggplot2)

all.slopes.mouse <- res.top.m$Segment.Slopes
colnames(all.slopes.mouse) <- paste0("mouse_slope", seq_len(ncol(all.slopes.mouse)))

all.bp.mouse <- res.top.m$Breakpoints
colnames(all.bp.mouse) <- paste0("mouse_breakpoint", seq_len(ncol(all.bp.mouse)))

all.slopes.human <- res.top.h$Segment.Slopes
colnames(all.slopes.human) <- paste0("human_slope", seq_len(ncol(all.slopes.human)))

all.bp.human <- res.top.h$Breakpoints
colnames(all.bp.human) <- paste0("human_breakpoint", seq_len(ncol(all.bp.human)))


all.bp.human <- (all.bp.human[ortho.genes.use$hgnc_symbol,])
all.slopes.human <- (all.slopes.human[ortho.genes.use$hgnc_symbol,])
all.bp.mouse <- (all.bp.mouse[ortho.genes.use$mgi_symbol,])
all.slopes.mouse <- (all.slopes.mouse[ortho.genes.use$mgi_symbol,])

X <- (c(all.bp.human[,1]))
Y <- (c(all.bp.mouse[,1]))

sum(X <= 250, na.rm=T) / length(X[!is.na(X)])
sum(Y <= 250, na.rm=T) / length(Y[!is.na(Y)])

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE, paired = TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)


pdf("PLOTS/histogram_firstBreakpoint_Orthologs_shortMouse-SuppFig1.pdf", height=2, width=3)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,300), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
dev.off()

library("yarrr")

X <- data.frame( Minute = X, Species = "Human")
Y <- data.frame( Minute = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_firstBreakpoint_Orthologs_shortMouse-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(1,2.5,1,.1), mgp=c(1.6,.5,0))
pirateplot(formula = Minute ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()


########################################################################################################
########################################################################################################
########################################################################################################

###################################################################################################
###################################################################################################
###################################################################################################

## Joint plots
########################################################################################################
## All GENES
# Supp Figure A and B and C. On Orthologs!
# Dynamic genes in common:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

sum(res.top.m$Breakpoint <= 60, na.rm=T)
sum(res.top.m$Breakpoint <= 100, na.rm=T)
sum(res.top.m$Breakpoint <= 250, na.rm=T)

sum(res.top.h$Breakpoint <= 60, na.rm=T)
sum(res.top.h$Breakpoint <= 100, na.rm=T)
sum(res.top.h$Breakpoint <= 250, na.rm=T)


X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library("yarrr")

X <- data.frame( Breakpoint = X, Species = "Human")
Y <- data.frame( Breakpoint = Y, Species = "Mouse")
longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_AllBreakpoints_anyGenes_shortMouse-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(1,2.5,1,.1), mgp=c(1.6,.5,0))
pirateplot(formula = Breakpoint ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()





# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.m <- timeUpDown.m[!is.na(timeUpDown.m)]
timeUpDown.h <- timeUpDown.h[!is.na(timeUpDown.h)]
sum(res.top.h$Trends[,1] !=0) / length(timeUpDown.h) # 92%
sum(res.top.m$Trends[,1] !=0) / length(timeUpDown.m) # 94%

propTestVals <- prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100

pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_mouseShort-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(mean(timeUpDown.h == 0), mean(timeUpDown.m==0), mean(timeUpDown.h > 0), mean(timeUpDown.m > 0))*100,
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="% Genes", xlab="",
ylim=c(0,100), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()



 
onlyDown.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == 1)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

propTestVals <- prop.test(c((all.mono.h), (all.mono.m)), c(nrow(res.top.h$Trends), nrow(res.top.m$Trends)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/numGenes_Monotonic_allTrendyGenes_mouseShort-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,2550), cex.axis=.6, cex.lab=.7)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()

all.mono.h / all.mono.m



########################################################################################################
########################################################################################################
########################################################################################################






# Ortholog Peaks:
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)

dim(ortho.peaks)

#### Step 1: Comparison of peak time in orthologs:
X <- c(subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$BreakPoint1)
Y <- c(subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$BreakPoint1)

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE, paired = TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library(ggplot2)
pdf("PLOTS/histogram_peakTime_orthologPeaks_mouseShort-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,30), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(300), cex=.6)
dev.off()


### Get peak times for all gene peaks now:

#### Step 1: Comparison of peak time in orthologs:
X <- c((peak_genes.h)$BreakPoint1)
Y <- c((peak_genes.m)$BreakPoint1)

length(X)
length(Y)

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library(ggplot2)
pdf("PLOTS/histogram_peakTime_anyPeaks_mouseShort-SuppFig1.pdf", height=2, width=2.5)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,200), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topleft', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(300), cex=.6)
dev.off()




###################################################################################################
## Slopes

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

library("yarrr")

# X <- data.frame( Slope = whichSlope.up.h, Species = "Human")
# Y <- data.frame( Slope = whichSlope.up.m, Species = "Mouse")
#
# longdata <- rbind(Y, X)
X <- data.frame( Slope = (whichSlope.up.h), Species = "Human")
Y <- data.frame( Slope = (whichSlope.up.m), Species = "Mouse")
longdata <- rbind(Y, X)
pdf("PLOTS/boxPlot_allPeak_slopeUP_mouseShort-SuppFig1.pdf", height=2, width=2.5)
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
pdf("PLOTS/boxPlot_allPeak_slopeDown_mouseShort-SuppFig1.pdf", height=2, width=2.5)
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
pdf("PLOTS/boxPlot_commonPeak_slopeUp_mouseShort-SuppFig1.pdf", height=2, width=2.5)
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
pdf("PLOTS/boxPlot_commonPeak_slopeDown_mouseShort-SuppFig1.pdf", height=2, width=2.5)
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
pdf("PLOTS/boxPlot_commonGenes_slopeUp_mouseShort-SuppFig1.pdf", height=2, width=2)
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
pdf("PLOTS/boxPlot_commonGenes_slopeDown_mouseShort-SuppFig1.pdf", height=2, width=2)
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()

