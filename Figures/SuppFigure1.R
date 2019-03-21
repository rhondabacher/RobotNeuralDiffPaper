setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")


## Get Common Peak

library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)


## Get Common UP

library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

upgenes.m <- c(names(which(res.top.m$Segment.Trends[,1] == 0 & res.top.m$Segment.Trends[,2] == 1)), 
                names(which(res.top.m$Segment.Trends[,1] == 1)))

upgenes.h <- c(names(which(res.top.h$Segment.Trends[,1] == 0 & res.top.h$Segment.Trends[,2] == 1)), 
                names(which(res.top.h$Segment.Trends[,1] == 1)))

set1 <- subset(ortho.genes.use, hgnc_symbol %in% upgenes.h)
ortho.upgenes <- subset(set1, mgi_symbol %in% upgenes.m)



#################################################################################################################
#################################################################################################################
#################################################################################################################

## Same as what is in Trendy package but wanted to control par directly. 
## Will add par control directly to the package soon.

fancyPlot2 <- function(DATA, tVectIn, trendyOutData, featureNames) {
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
	   }}else {
							   IDseg <- ID[1]
						       useCol <- switch(names(which.max(table(IDseg))), 
						       "0" = "black", 
						       "-1" = "cornflowerblue", 
						       "1" = "coral1")
							   	lines(tVectIn, FIT, lwd = 5, col=useCol)
						   }
						 
		 
	   par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
	       mar = c(0, 0, 4, 0), new = TRUE)
	   plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

	
}
################################################################################################################
################################################################################################################

# Up trend examples:

pdf("PLOTS/Mouse_OrthologGenes_ExpressionScatter_HSPA8_SuppFig1.pdf", height=4.5, width=6.5, useDingbats=FALSE)
par(mfrow=c(1,1), cex=1.5, cex.lab=1.2, cex.axis=1.2, cex.main=1.1, mar=c(4,4,2,1), mgp=c(2.5,1,0))
XX<- fancyPlot2(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Hspa8", 
            trendyOutData = seg.mouse)
dev.off()

pdf("PLOTS/Mouse_OrthologGenes_ExpressionScatter_HSPA14_uppFig1.pdf", height=4.5, width=6.5, useDingbats=FALSE)
par(mfrow=c(1,1), cex=1.5, cex.lab=1.2, cex.axis=1.2, cex.main=1.1, mar=c(4,4,2,1), mgp=c(2.5,1,0))
XX<- fancyPlot2(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Hspa14", 
            trendyOutData = seg.mouse)
dev.off()


pdf("PLOTS/Mouse_OrthologGenes_ExpressionScatter_IVNS1ABP_SuppFig1.pdf", height=4.5, width=6.5, useDingbats=FALSE)
par(mfrow=c(1,1), cex=1.5, cex.lab=1.2, cex.axis=1.2, cex.main=1.1, mar=c(4,4,2,1), mgp=c(2.5,1,0))
XX<- fancyPlot2(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Ivns1abp", 
            trendyOutData = seg.mouse)
dev.off()


pdf("PLOTS/Mouse_OrthologGenes_ExpressionScatter_MCM7_SuppFig1.pdf", height=4.5, width=6.5, useDingbats=FALSE)
par(mfrow=c(1,1), cex=1.5, cex.lab=1.2, cex.axis=1.2, cex.main=1.1, mar=c(4,4,2,1), mgp=c(2.5,1,0))
XX<- fancyPlot2(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Mcm7", 
            trendyOutData = seg.mouse)
dev.off()



# Time of max expression for Up genes:
peak.com.m <- ortho.upgenes$mgi_symbol
whichTime.up.m<-c()
for(i in 1:length(peak.com.m)) {
	keepT <- which(res.top.m$Segment.Trends[peak.com.m[i],] != 1)[1]
	
	if(is.na(keepT)) {
		whichTime.up.m[i] <- 600
	} else if (keepT==1) {
		keepT <- which(res.top.m$Segment.Trends[peak.com.m[i],] == 1)[1]
		whichTime.up.m[i] <- res.top.m$Breakpoints[peak.com.m[i],(keepT)]
		if (is.na(whichTime.up.m[i])) {whichTime.up.m[i] <- 600}
	} else {
	whichTime.up.m[i] <- res.top.m$Breakpoints[peak.com.m[i],(keepT-1)]
	}
}
names(whichTime.up.m) <- peak.com.m


peak.com.h <- set1$hgnc_symbol
whichTime.up.h<-c()
for(i in 1:length(peak.com.h)) {
	keepT <- which(res.top.h$Segment.Trends[peak.com.h[i],] != 1)[1]
	
	if(is.na(keepT)) {
		whichTime.up.h[i] <- 600
	} else if (keepT==1) {
		keepT <- which(res.top.h$Segment.Trends[peak.com.h[i],] == 1)[1]
		whichTime.up.h[i] <- res.top.h$Breakpoints[peak.com.h[i],(keepT)]
		if (is.na(whichTime.up.h[i])) {whichTime.up.h[i] <- 600}
	} else {
	whichTime.up.h[i] <- res.top.h$Breakpoints[peak.com.h[i],(keepT-1)]
	}
}
names(whichTime.up.h) <- peak.com.h


X = whichTime.up.h
Y = whichTime.up.m



pdf("PLOTS/spectrum_EndTimeOfUp_Orthologs.pdf", height=5, width=18)
par(mar=c(6,3,1,1), mfrow=c(1,1), mgp=c(1,4,0))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', ylab="", xlab="",main="", cex.axis=5)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=5, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=5, col="cornflowerblue")
dev.off()



