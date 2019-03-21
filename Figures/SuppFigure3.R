setwd("~/RobotSeq/")
library(Trendy)



load("RDATA/trendy_run_mouse_InVitro.Rdata")
seg.mouse <- ready.trendy.scaled
data.norm.m <- data.norm
data.norm.scale.m <- data.norm.scale
t.v.m <- t.v

load("OUT2/RDATA/run_trendy_inVitro_mouse_both.Rdata")
seg.human <- ready.trendy.scaled
data.norm.h <- data.norm
data.norm.scale.h <- data.norm.scale
t.v.h <- t.v


# First get top genes in both species, per Trendy Fig 2 use higher adj. R^2 for shorter time course:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .8)
res.top.h <- topTrendy(seg.human, .8)

trendy.summary.m <- formatResults(res.top.m)
trendy.summary.h <- formatResults(res.top.h)

peak_genes.h <- extractPattern(seg.human, Pattern = c("up","down"), adjR2Cut=.8)
peak_genes.m <- extractPattern(seg.mouse, Pattern = c("up","down"), adjR2Cut=.8)

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

top.mouse <- top.mouse[peak_genes.m[,1],]
top.human <- top.human[peak_genes.h[,1],]

orth.genes.clean1 <- orth.genes[which(orth.genes[,1]!="" & orth.genes[,2]!=""),]
orth.genes.clean <- orth.genes.clean1[!duplicated(orth.genes.clean1),]
top1 <- merge(orth.genes.clean, top.human, by="hgnc_symbol")
top2 <- merge(top1,top.mouse, by="mgi_symbol")
top2 <- top2[!duplicated(top2),]
dupg <-top2[which(duplicated(top2[,2])),2] 
subset(top2, hgnc_symbol %in% dupg)
dupg <-top2[which(duplicated(top2[,1])),2]
subset(top2, hgnc_symbol %in% dupg)

# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse, !(Gene %in% top2$mgi_symbol))
toAdd <- intersect(toupper(mouse.genes.check1[,1]), top.human$Gene)
toAdd <- subset(top.mouse, toupper(Gene) %in% toAdd)[,1]

recovered.g <- data.frame(mgi_symbol = c(toAdd), hgnc_symbol = toupper(toAdd), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)
dim(top2)

# Get peak up slope and down slope

peak.com.h <- subset(peak_genes.h, Gene %in% top2$hgnc_symbol)
whichSlope.up.h<-c()
whichSlope.down.h<-c()
for(i in 1:nrow(peak.com.h)) {
	keep <- which(round(res.top.h$Breakpoints[peak.com.h[i,1],]) == round(c(peak.com.h[i,2])))
	whichSlope.up.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],keep]
	whichSlope.down.h[i] <- res.top.h$Segment.Slopes[peak.com.h[i,1],(keep+1)]
	
}
names(whichSlope.up.h) <- peak.com.h[,1]
names(whichSlope.down.h) <- peak.com.h[,1]

peak.com.m <- subset(peak_genes.m, Gene %in% top2$mgi_symbol)
whichSlope.up.m<-c()
whichSlope.down.m<-c()
for(i in 1:nrow(peak.com.m)) {
	keep <- which(round(res.top.m$Breakpoints[peak.com.m[i,1],]) == round(c(peak.com.m[i,2])))
	whichSlope.up.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],keep]
	whichSlope.down.m[i] <- res.top.m$Segment.Slopes[peak.com.m[i,1],(keep+1)]	
}
names(whichSlope.up.m) <- peak.com.m[,1]
names(whichSlope.down.m) <- peak.com.m[,1]


X <- whichSlope.up.h
Y <- whichSlope.up.m

wilcox.test(X,Y)
PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
  PP
library(ggplot2)
pdf("PLOTS/histogram_commonPeak_slopeUP_INVITRO_SuppFigure3.pdf", height=6, width=7.5)
par(mar=c(6,6,3,1))
hist(X, xlim=c(0,.55), ylim=c(0,15), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, .55, length.out=50), 
	main="", xlab="Slope (Scaled Expression/Minute)",
	cex.axis=1.8, cex.lab=2)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, .55, length.out=50))
legend('topright', c("Mouse","Human"), lwd=3, col=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2)
dev.off()


library("yarrr")

X <- data.frame( Slope = whichSlope.up.h, Species = "Human")
Y <- data.frame( Slope = whichSlope.up.m, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_commonPeak_slopeUP_INVITRO_SuppFigure3.pdf", height=6, width=6)
par(mar=c(5,7,2,1), mgp = c(5, .5, 0))
pirateplot(formula = Slope ~ Species,
           data = longdata,
           xlab = "",
           ylab = "Slope", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
dev.off()



X <- whichSlope.down.h
Y <- whichSlope.down.m

wilcox.test(X,Y)
PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
  PP
library(ggplot2)
pdf("PLOTS/histogram_commonPeak_slopeDown_INVITRO_SuppFigure3.pdf", height=6, width=7.5)
par(mar=c(6,6,3,1))
hist(X, xlim=c(-.45,0), ylim=c(0,20), border="brown3",
	col=alpha("brown1", .6), breaks = seq(-.45, 0, length.out=50), 
	main="", xlab="Slope (Scaled Expression/Minute)",
	cex.axis=1.8, cex.lab=2)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(-.45, 0, length.out=50))
legend('topleft', c("Mouse","Human"), lwd=3, col=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2)
dev.off()


library("yarrr")

X <- data.frame( Slope = whichSlope.down.h, Species = "Human")
Y <- data.frame( Slope = whichSlope.down.m, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_commonPeak_slopeDown_INVITRO_SuppFigure3.pdf", height=6, width=6)
par(mar=c(5,8,2,1), mgp = c(6, .5, 0))
pirateplot(formula = Slope ~ Species,
           data = longdata,
           xlab = "",
           ylab = "Slope", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
dev.off()



