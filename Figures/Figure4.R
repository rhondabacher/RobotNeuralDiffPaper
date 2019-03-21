setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

library(Trendy)

# Ortholog genes:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

## Tables for supplement:

tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[names(sort(res.top.m$Breakpoints[ortho.genes.use$mgi_symbol,1],na.last=TRUE)), ]
write.table(tosave.m, file="TABLES/summary_mouse_allOrthologs.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[names(sort(res.top.h$Breakpoints[ortho.genes.use$hgnc_symbol,1],na.last=TRUE)), ]
write.table(tosave.h, file="TABLES/summary_human_allOrthologs.csv", row.names=F, quote=F, sep=",")


########################################################################################################
## Get breakpoints and slopes for ortholog genes

library(gplots)

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

X <- (na.omit(c(all.bp.human[,1])))
Y <- (na.omit(c(all.bp.mouse[,1])))

PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
  PP
  sum(X <= 250) / length(X)
  sum(Y <= 250) / length(Y)
  
pdf("PLOTS/histogram_firstBreakpoint_Orthologs_Figure4.pdf", height=8, width=12)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,15), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=2.5, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2)
dev.off()

pdf("PLOTS/spectrum_firstBreakpoint_Orthologs_Figure4.pdf", height=3.5, width=16)
par(mar=c(3,3,1,1), mfrow=c(1,1), mgp=c(0,2,0))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', ylab="", xlab="",main="", cex.axis=3)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=2, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=2, col="cornflowerblue")
dev.off()

library("yarrr")

X <- data.frame( Minute = X, Species = "Human")
Y <- data.frame( Minute = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_firstBreakpoint_Orthologs_Figure4.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
pirateplot(formula = Minute ~ Species,
           data = longdata,
           xlab = "",
           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()


########################################################################################################

pdf("PLOTS/allPossibleOrthologs_Human.pdf", height=6, width=6)
par(mfrow=c(2,2))
plotFeature(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = "IFFO2", legendLocation='bottom',
            trendyOutData = seg.human)
dev.off()

pdf("PLOTS/allPossibleOrthologs_Mouse.pdf", height=6, width=6)
par(mfrow=c(2,2))
plotFeature(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Iffo2", legendLocation='bottom',
            trendyOutData = seg.mouse)
dev.off()

#######################################################################################################################
##########################################################################################################################

## Make the plots for peak Orthologs:
								
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)


## Tables for supplement:
tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$Gene, ]
write.table(tosave.m, file="TABLES/summary_mouse_peakOrthologs.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$Gene, ]
write.table(tosave.h, file="TABLES/summary_human_peakOrthologs.csv", row.names=F, quote=F, sep=",")



## Only for peak orthologs!
X <- c(subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$BreakPoint1)
Y <- c(subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$BreakPoint1)

sum(X <= 180) / length(X)
sum(Y <= 180) / length(Y)

sum(X > 300) / length(X)
sum(Y > 300) / length(Y)


PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
  PP
pdf("PLOTS/histogram_firstBreakpoint_orthologPeaks_Figure4.pdf", height=8, width=12)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,15), border="brown3", 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=2.5, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6),border="dodgerblue3",  breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6),alpha("brown1", .6)), cex=2)
dev.off()

pdf("PLOTS/spectrum_firstBreakpoint_orthologPeaks_Figure4.pdf", height=3.5, width=16)
par(mar=c(3,3,1,1), mfrow=c(1,1), mgp=c(0,2,0))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', ylab="", xlab="",main="", cex.axis=3)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=2, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=2, col="cornflowerblue")
dev.off()

library("yarrr")

X <- data.frame( Minute = X, Species = "Human")
Y <- data.frame( Minute = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_firstBreakpoint_orthologsPeaks_Figure4.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
pirateplot(formula = Minute ~ Species,
           data = longdata,
           xlab = "",
           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()





################################################################################################################
################################################################################################################

# Now make some individual gene plots:

##################################################################################################################
##################################################################################################################
## Same as what is in Trendy package but wanted to control par directly. 
## Will add par control directly to the package soon.
fancyPlot <- function(DATA, tVectIn, trendyOutData, featureNames) {
	plot(tVectIn, DATA[featureNames,], pch=20, col="#696969", main=featureNames, ylab="Scaled Expression", xlab="Minute", 
	           cex.axis=1.4, cex.lab=1.5, xaxt='n', yaxt='n', cex.main=2)
  axis(1, at=c(0,200,400,600), cex.lab=1.5, cex.axis=1.5)
  axis(2, at=c(0,.5,1), cex.lab=1.5, cex.axis=1.5)
	trendyOutData <- trendyOutData[[featureNames]]
	lines(tVectIn, trendyOutData$Fitted.Values, lwd = 3, col="#ededed")
	abline(v = trendyOutData$Breakpoints, lty = 2, lwd = 3, col="chartreuse3")
	ID <- trendyOutData$Trends
	FIT <- trendyOutData$Fitted.Values
	BKS <- c(0, trendyOutData$Breakpoints, max(tVectIn))
	if (length(BKS) > 3 | (length(BKS) == 3 & !is.na(BKS[2]))) {
	   for (i in seq_len(length(BKS) - 1)) {
	       toCol <- which(tVectIn <= BKS[i+1] & tVectIn >= BKS[i])
	       IDseg <- ID[toCol]
	       useCol <- switch(names(which.max(table(IDseg))), 
	       "0" = "black", 
	       "-1" = "cornflowerblue", 
	       "1" = "coral1")
	       lines(tVectIn[toCol], FIT[toCol], lwd = 5, col=useCol)
	   }


	}
}

################################################################################
################################################################################
################################################################################

m1 <- extractPattern(seg.mouse, .5, Pattern = c("up", "down"), 0)
h1 <- extractPattern(seg.human, .5, Pattern = c("up", "down"), 0)
X1 <- ortho.peaks[which(ortho.peaks[,2] %in% h1$Gene & ortho.peaks[,1] %in% m1$Gene),]

## Plot first 10
toplot1 <- subset(m1, Gene %in% X1[,1])
toplot1 <- toplot1[1:10,1]
toplot1.h <- toupper(toplot1)

for(i in 1:10){

pdf(paste0("PLOTS/OrthologGenes_ExpressionScatter_",toplot1.h[i],".pdf"), height=9, width=6.5, useDingbats=FALSE)
par(mfrow=c(2,1), cex=1.5, cex.lab=1.2, cex.axis=1.2, cex.main=1.1, mar=c(4,4,2,1), mgp=c(2.5,1,0))
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
dev.off()

}



