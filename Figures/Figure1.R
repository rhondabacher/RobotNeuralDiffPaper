setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")


######################################################################################################
######################################################################################################

## Same as what is in Trendy package but wanted to control par directly. 
## Will add par control directly to the package soon.

fancyPlot <- function(DATA, tVectIn, trendyOutData, featureNames) {
	plot(tVectIn, DATA[featureNames,], pch=20, col="#696969", cex=.6,
  main=featureNames, ylab="Scaled Expression", xlab="Minute", 
	            cex.lab=1, xaxt='n', yaxt='n', cex.main=1)
  axis(1, at=c(0,200,400,600), cex.axis=1)
  axis(2, at=c(0,.5,1), cex.axis=1)
  if(!is.null(trendyOutData)) {
	trendyOutData <- trendyOutData[[featureNames]]
	lines(tVectIn, trendyOutData$Fitted.Values, lwd = 1, col="#ededed")
	abline(v = trendyOutData$Breakpoints, lty = 2, lwd = 1, col="chartreuse3")
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
	       lines(tVectIn[toCol], FIT[toCol], lwd = 1, col=useCol)
	   }} else {
							   IDseg <- ID[1]
						       useCol <- switch(names(which.max(table(IDseg))), 
						       "0" = "black", 
						       "-1" = "cornflowerblue", 
						       "1" = "coral1")
							   	lines(tVectIn, FIT, lwd = 1, col=useCol)
						   }
             }

}
################################################################################################################ 
################################################################################################################
pdf("PLOTS/Figure1_GeneScatter.pdf", height=3, width=7.5, useDingbats=F)
par(mfrow=c(2,6), mar=c(3,2,2,.1), mgp=c(1.1,.4,0))

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Sept7", 
            trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Apex1", 
            trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Nid1", 
            trendyOutData = seg.mouse)
 
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
                        featureNames = "Lefty1", 
                        trendyOutData = seg.mouse)  

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Eml1", 
            trendyOutData = seg.mouse)     
						            
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Tpm1", 
            trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
                        featureNames = "Exo1", 
                        trendyOutData = seg.mouse)  

fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
          featureNames = "Myc", 
          trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Gpc4", 
            trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
              featureNames = "Bclaf1", 
              trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Slc30a1", 
            trendyOutData = seg.mouse)
fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = "Fgfbp1", 
            trendyOutData = seg.mouse)
dev.off()


#####################################################################
#####################################################################

# How many patterns in each species for Figure 1B?
library(Trendy)

## MOUSE FIRST:
all.top <- topTrendy(seg.mouse, .2)
print(length(all.top$A))

# Not Dynamic:
length(seg.mouse) - length(all.top$A)

all.total <- topTrendy(seg.mouse, -1)
trendDirs <- apply(all.total$Segment.Trends, 1, function(x) {
  length(unique(x[!is.na(x)]))
})
firstDir <- all.total$Segment.Trends[,1]
sameGenes <- (names(which(trendDirs == 1 & firstDir == 0)))
length(sameGenes)


# Monotonic
trendDirs <- apply(all.top$Segment.Trends, 1, function(x) {
  length(unique(x[!is.na(x)]))
})
firstDir <- all.top$Segment.Trends[,1]
upGenes <- names(which(trendDirs == 1 & firstDir == 1))
downGenes <- names(which(trendDirs == 1 & firstDir == -1))

length(upGenes)
length(downGenes)

# Delay Peak/Dip
delayPeakGenes <- (extractPattern(seg.mouse, Pattern=c("same", "up", "down"), adjR2Cut =.2))
# delaySlowPeakGenes <- (extractPattern(seg.mouse, Pattern=c("same", "up", "same", "down"), adjR2Cut =.2))

delayDipGenes <- (extractPattern(seg.mouse, Pattern=c("same", "down", "up"), adjR2Cut =.2))
delaySlowDipGenes <- (extractPattern(seg.mouse, Pattern=c("same", "down", "same", "up"), adjR2Cut =.2))
delayDipGenes <- rbind(delayDipGenes, delaySlowDipGenes[,1:3])

delayPeakGenes <- subset(delayPeakGenes, delayPeakGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0)))[,1]
delayDipGenes <- subset(delayDipGenes, delayDipGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0)))[,1]

length(delayPeakGenes) ## 
length(delayDipGenes) ## 

# Cyclic in either way:
cycleGenes1 <- (extractPattern(seg.mouse, Pattern=c("up","down","up"), adjR2Cut =.2))
cycleGenes2 <- (extractPattern(seg.mouse, Pattern=c("up","same","down","up"), adjR2Cut =.2))
cycleGenes3 <- (extractPattern(seg.mouse, Pattern=c("up","down","same","up"), adjR2Cut =.2))
cycleGenes1 <- rbind(cycleGenes1, cycleGenes2[,1:3], cycleGenes3[,1:3])
cycleGenes2 <- (extractPattern(seg.mouse, Pattern=c("down","up","down"), adjR2Cut =.2))
cycleGenes3 <- (extractPattern(seg.mouse, Pattern=c("down","same","up","down"), adjR2Cut =.2))
cycleGenes4 <- (extractPattern(seg.mouse, Pattern=c("down","up","same","down"), adjR2Cut =.2))
cycleGenes1 <- rbind(cycleGenes1, cycleGenes2[,1:3], cycleGenes3[,1:3], cycleGenes4[,1:3])

cycleGenes5 <- (extractPattern(seg.mouse, Pattern=c("up","down","down","up"), adjR2Cut =.2))
cycleGenes6 <- (extractPattern(seg.mouse, Pattern=c("down","up","up","down"), adjR2Cut =.2))
cycleGenes7 <- (extractPattern(seg.mouse, Pattern=c("up","down","same","down","up"), adjR2Cut =.2))
cycleGenes1 <- rbind(cycleGenes1, cycleGenes5[,1:3], cycleGenes6[,1:3], cycleGenes7[,1:3])
cycleGenes <- unique(setdiff(cycleGenes1[,1], c(delayPeakGenes, delayDipGenes)))
length(cycleGenes)

# Peak only:
peakGenes <- extractPattern(seg.mouse, Pattern=c("up","down"), adjR2Cut =.2)
slowPeakGenes <- extractPattern(seg.mouse, Pattern=c("up","same","down"), adjR2Cut =.2)
onlyPeak <- unique(setdiff(c(peakGenes[,1], slowPeakGenes[,1]), c(cycleGenes, delayPeakGenes, delayDipGenes)))
length(onlyPeak)

dipGenes <- extractPattern(seg.mouse, Pattern=c("down", "up"), adjR2Cut =.2)
slowDipGenes <- extractPattern(seg.mouse, Pattern=c("down","same","up"), adjR2Cut =.2) 
allDip <- unique(setdiff(c(dipGenes[,1], slowDipGenes[,1]), 
    c(cycleGenes, delayPeakGenes, delayDipGenes)))
length(allDip)

# Delay Up/Down
delayupGenes <- extractPattern(seg.mouse, Pattern=c("same", "up"), adjR2Cut =.2)
delayupGenes1 <- subset(delayupGenes, delayupGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & all.top$Segment.Trends[,2] == 1)))
delayupGenes2 <- subset(delayupGenes, delayupGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & 
                                                                  all.top$Segment.Trends[,2] == 0 & all.top$Segment.Trends[,3] == 1)))
delayupGenes <- rbind(delayupGenes1, delayupGenes2)
delayUp <- setdiff(delayupGenes[,1], c(allDip, onlyPeak, delayPeakGenes, delayDipGenes))
length(delayUp)

delaydownGenes <- extractPattern(seg.mouse, Pattern=c("same", "down"), adjR2Cut =.2)
delaydownGenes1 <- subset(delaydownGenes, delaydownGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & all.top$Segment.Trends[,2] == -1)))
delaydownGenes2 <- subset(delaydownGenes, delaydownGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & 
                                                                  all.top$Segment.Trends[,2] == 0 & all.top$Segment.Trends[,3] == -1)))
delaydownGenes <- rbind(delaydownGenes1, delaydownGenes2)
delayDown <- setdiff(delaydownGenes[,1], c(allDip, onlyPeak, delayPeakGenes, delayDipGenes))
length(delayDown)

alldynamic <- c(delayUp, delayDown, allDip, 
onlyPeak, cycleGenes, delayDipGenes, delayPeakGenes, downGenes, upGenes)

check1 <- setdiff(names(all.top$AdjustedR2), alldynamic)
length(check1)

immediateOn <- extractPattern(seg.mouse, Pattern=c("up", "same"), adjR2Cut =.2)
length(setdiff(immediateOn[,1], alldynamic))

immediateDown <- extractPattern(seg.mouse, Pattern=c("down", "same"), adjR2Cut =.2)
length(setdiff(immediateDown[,1], alldynamic))

other <- base::setdiff(names(all.top$AdjustedR2), c(alldynamic, immediateDown[,1], immediateOn[,1]))
length(other)


length(alldynamic)
length(unique(alldynamic))

alldynamic[which(duplicated(alldynamic))]




tosave <- data.frame(Pattern = c("TotalNotDynamic", "TotalDynamic", "MonotonicUp", 
																	"MonotonicDown", "ImmediateOff", "ImmediateOn",
																	"Peak", "Dip", "Cyclic", "DelayUp", "DelayDown", 
																	"DelayPeak", "DelayDip", "Unclassified"),
										 NumGenes = c(length(seg.mouse) - length(all.top$A), 
										 							length(all.top$A), length(upGenes), length(downGenes), 
																	length(immediateOn[,1]), length(immediateDown[,1]),
																	length(onlyPeak), length(allDip), length(cycleGenes), 
																	length(delayUp), length(allDip), 
																	length(delayPeakGenes), length(delayDipGenes), 
																	length(other))
)

write.table(tosave, file="TABLES/NumberOfGenes_PerDyanmicPattern_Mouse.csv", row.names=F, quote=F, sep=",")


library(rowr)
tosave <- cbind.fill(
  MonotonicUp = upGenes, 
  MonotonicDown =  downGenes, 
  ImmediateOff = immediateOn[,1], 
  ImmediateOn = immediateDown[,1],
  Peak = onlyPeak, 
  Dip = allDip, 
  Cyclic = cycleGenes, 
  DelayUp = delayUp, 
  DelayDown = allDip, 
  DelayPeak = delayPeakGenes,
  DelayDip = delayDipGenes, 
  Unclassified = other,
  fill="")
colnames(tosave) <- c("MonotonicUp", "MonotonicDown",
  "ImmediateOff", "ImmediateOn",
  "Peak", "Dip", "Cyclic",
  "DelayUp", "DelayDown", "DelayPeak",
  "DelayDip", "Unclassified")
str(tosave)
write.table(tosave, file="TABLES/ListOfGenes_PerDyanmicPattern_Mouse.csv", row.names=F, quote=F, sep=",")







rm(list=ls())



load("RDATA/jointPlots_loadDataBoth.Rdata")


## HUMAN
library(Trendy)

## MOUSE FIRST:
all.top <- topTrendy(seg.human, .2)
print(length(all.top$A))

# Not Dynamic:
length(seg.human) - length(all.top$A)

all.total <- topTrendy(seg.human, -1)
trendDirs <- apply(all.total$Segment.Trends, 1, function(x) {
  length(unique(x[!is.na(x)]))
})
firstDir <- all.total$Segment.Trends[,1]
sameGenes <- (names(which(trendDirs == 1 & firstDir == 0)))
length(sameGenes)


# Monotonic
trendDirs <- apply(all.top$Segment.Trends, 1, function(x) {
  length(unique(x[!is.na(x)]))
})
firstDir <- all.top$Segment.Trends[,1]
upGenes <- names(which(trendDirs == 1 & firstDir == 1))
downGenes <- names(which(trendDirs == 1 & firstDir == -1))

length(upGenes)
length(downGenes)

# Delay Peak/Dip
delayPeakGenes <- (extractPattern(seg.human, Pattern=c("same", "up", "down"), adjR2Cut =.2))
delaySlowPeakGenes <- (extractPattern(seg.human, Pattern=c("same", "up", "same", "down"), adjR2Cut =.2))

delayDipGenes <- (extractPattern(seg.human, Pattern=c("same", "down", "up"), adjR2Cut =.2))
delaySlowDipGenes <- (extractPattern(seg.human, Pattern=c("same", "down", "same", "up"), adjR2Cut =.2))
delayDipGenes <- rbind(delayDipGenes, delaySlowDipGenes[,1:3])

delayPeakGenes <- subset(delayPeakGenes, delayPeakGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0)))[,1]
delayDipGenes <- subset(delayDipGenes, delayDipGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0)))[,1]

length(delayPeakGenes) ## 
length(delayDipGenes) ## 

# Cyclic in either way:
cycleGenes1 <- (extractPattern(seg.human, Pattern=c("up","down","up"), adjR2Cut =.2))
cycleGenes2 <- (extractPattern(seg.human, Pattern=c("up","same","down","up"), adjR2Cut =.2))
cycleGenes3 <- (extractPattern(seg.human, Pattern=c("down","up","same","up", "down"), adjR2Cut =.2))
cycleGenes1 <- rbind(cycleGenes1, cycleGenes2[,1:3], cycleGenes3[,1:3])
cycleGenes2 <- (extractPattern(seg.human, Pattern=c("down","up","down"), adjR2Cut =.2))
cycleGenes3 <- (extractPattern(seg.human, Pattern=c("down","same","up","down"), adjR2Cut =.2))
cycleGenes4 <- (extractPattern(seg.human, Pattern=c("down","up","same","down"), adjR2Cut =.2))
cycleGenes5 <- (extractPattern(seg.human, Pattern=c("down","up","up","down"), adjR2Cut =.2))
cycleGenes6 <- (extractPattern(seg.human, Pattern=c("up","down","down","up"), adjR2Cut =.2))
cycleGenes1 <- rbind(cycleGenes1, cycleGenes2[,1:3], cycleGenes3[,1:3], cycleGenes4[,1:3]
  , cycleGenes5[,1:3], cycleGenes6[,1:3])
cycleGenes <- unique(setdiff(cycleGenes1[,1], c(delayPeakGenes, delayDipGenes)))
length(cycleGenes)

# Peak only:
peakGenes <- extractPattern(seg.human, Pattern=c("up","down"), adjR2Cut =.2)
slowPeakGenes <- extractPattern(seg.human, Pattern=c("up","same","down"), adjR2Cut =.2)
onlyPeak <- unique(setdiff(c(peakGenes[,1], slowPeakGenes[,1]), c(cycleGenes, delayPeakGenes, delayDipGenes)))
length(onlyPeak)

dipGenes <- extractPattern(seg.human, Pattern=c("down", "up"), adjR2Cut =.2)
slowDipGenes <- extractPattern(seg.human, Pattern=c("down","same","up"), adjR2Cut =.2) 
allDip <- unique(setdiff(c(dipGenes[,1], slowDipGenes[,1]), 
    c(cycleGenes, delayPeakGenes, delayDipGenes)))
length(allDip)

# Delay Up/Down
delayupGenes <- extractPattern(seg.human, Pattern=c("same", "up"), adjR2Cut =.2)
delayupGenes1 <- subset(delayupGenes, delayupGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & all.top$Segment.Trends[,2] == 1)))
delayupGenes2 <- subset(delayupGenes, delayupGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & 
                                                                  all.top$Segment.Trends[,2] == 0 & all.top$Segment.Trends[,3] == 1)))
delayupGenes <- rbind(delayupGenes1, delayupGenes2)
delayUp <- setdiff(delayupGenes[,1], c(allDip, onlyPeak, delayPeakGenes, delayDipGenes))
length(delayUp)

delaydownGenes <- extractPattern(seg.human, Pattern=c("same", "down"), adjR2Cut =.2)
delaydownGenes1 <- subset(delaydownGenes, delaydownGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & all.top$Segment.Trends[,2] == -1)))
delaydownGenes2 <- subset(delaydownGenes, delaydownGenes[,1] %in% names(which(all.top$Segment.Trends[,1] == 0 & 
                                                                  all.top$Segment.Trends[,2] == 0 & all.top$Segment.Trends[,3] == -1)))
delaydownGenes <- rbind(delaydownGenes1, delaydownGenes2)
delayDown <- setdiff(delaydownGenes[,1], c(allDip, onlyPeak, delayPeakGenes, delayDipGenes))
length(delayDown)

alldynamic <- c(delayUp, delayDown, allDip, 
onlyPeak, cycleGenes, delayDipGenes, delayPeakGenes, downGenes, upGenes)

check1 <- setdiff(names(all.top$AdjustedR2), alldynamic)
length(check1)

immediateOn <- extractPattern(seg.human, Pattern=c("up", "same"), adjR2Cut =.2)
length(setdiff(immediateOn[,1], alldynamic))

immediateDown <- extractPattern(seg.human, Pattern=c("down", "same"), adjR2Cut =.2)
length(setdiff(immediateDown[,1], alldynamic))

other <- base::setdiff(names(all.top$AdjustedR2), c(alldynamic, immediateDown[,1], immediateOn[,1]))
length(other)


length(alldynamic)
length(unique(alldynamic))

alldynamic[which(duplicated(alldynamic))]




tosave <- data.frame(Pattern = c("TotalNotDynamic", "TotalDynamic", "MonotonicUp", 
																	"MonotonicDown", "ImmediateOff", "ImmediateOn",
																	"Peak", "Dip", "Cyclic", "DelayUp", "DelayDown", 
																	"DelayPeak", "DelayDip", "Unclassified"),
										 NumGenes = c(length(seg.human) - length(all.top$A), 
										 							length(all.top$A), length(upGenes), length(downGenes), 
																	length(immediateOn[,1]), length(immediateDown[,1]),
																	length(onlyPeak), length(allDip), length(cycleGenes), 
																	length(delayUp), length(allDip), 
																	length(delayPeakGenes), length(delayDipGenes), 
																	length(other))
)

write.table(tosave, file="TABLES/NumberOfGenes_PerDyanmicPattern_Human.csv", row.names=F, quote=F, sep=",")


library(rowr)
tosave <- cbind.fill(
  MonotonicUp = upGenes, 
  MonotonicDown =  downGenes, 
  ImmediateOff = immediateOn[,1], 
  ImmediateOn = immediateDown[,1],
  Peak = onlyPeak, 
  Dip = allDip, 
  Cyclic = cycleGenes, 
  DelayUp = delayUp, 
  DelayDown = allDip, 
  DelayPeak = delayPeakGenes,
  DelayDip = delayDipGenes, 
  Unclassified = other,
  fill="")
colnames(tosave) <- c("MonotonicUp", "MonotonicDown",
  "ImmediateOff", "ImmediateOn",
  "Peak", "Dip", "Cyclic",
  "DelayUp", "DelayDown", "DelayPeak",
  "DelayDip", "Unclassified")
str(tosave)
write.table(tosave, file="TABLES/ListOfGenes_PerDyanmicPattern_Human.csv", row.names=F, quote=F, sep=",")




