setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")



#################################################################################################################
#################################################################################################################

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
	   }} else {
							   IDseg <- ID[1]
						       useCol <- switch(names(which.max(table(IDseg))), 
						       "0" = "black", 
						       "-1" = "cornflowerblue", 
						       "1" = "coral1")
							   	lines(tVectIn, FIT, lwd = 5, col=useCol)
						   }

}
################################################################################################################
################################################################################################################

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Wsb1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Tpm1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Slc30a1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Eml1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Bclaf1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Nid1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Gpc4", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Apex1", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Sept7", 
            trendyOutData = seg.mouse)

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Myc", 
            trendyOutData = seg.mouse)


			