setwd("~/RobotNeuralDiffPaper/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

rm(list=setdiff(ls(), "orth.genes.clean"))

library(Trendy)
load("RDATA/trendy_run_mouse_InVitro.Rdata")
seg.mouse <- results(seg.all.scaled)
data.norm.m <- data.norm
data.norm.scale.m <- data.norm.scale
t.v.m <- t.v

load("RDATA/trendy_run_human_InVitro.Rdata")
seg.human <- results(seg.all.scaled)
data.norm.h <- data.norm
data.norm.scale.h <- data.norm.scale
t.v.h <- t.v



# First get top genes in both species, per Trendy Fig 2 use higher adj. R^2 for shorter time course:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

trendy.summary.m <- formatResults(res.top.m)
trendy.summary.h <- formatResults(res.top.h)

peak_genes <- extractPattern(seg.human, Pattern=c("up", "down"), adjR2Cut =.5)
peak_genes <- rbind(peak_genes, extractPattern(seg.human, Pattern=c("up", "same", "down"), adjR2Cut =.5)[,1:2])
peak_genes.h <- peak_genes

peak_genes <- extractPattern(seg.mouse, Pattern=c("up", "down"), adjR2Cut =.5)
peak_genes <- rbind(peak_genes, extractPattern(seg.mouse, Pattern=c("up", "same", "down"), adjR2Cut =.5)[,1:2])
peak_genes.m <- peak_genes

peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]


top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), 
                        mgi_symbol=names(res.top.m$AdjustedR2), 
                        row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), 
                        hgnc_symbol=names(res.top.h$AdjustedR2), 
                        row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)


##############################################################################
##############################################################################
# Now compare the up and down slope for genes that are otholog peaks:
##############################################################################


peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]



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

# Get peak up slope and down slope



set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)


peak.com.h <- subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)
whichSlope.up.h <- c()
whichSlope.down.h <- c()
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
library(yarrr)
pdf("PLOTS/boxPlot_commonPeak_slopeUp_Barry2017.pdf", height=2, width=2)
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
pdf("PLOTS/boxPlot_commonPeak_slopeDown_Barry2017.pdf", height=2, width=2)
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()



