setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)


# Ortholog Peaks:
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)

dim(ortho.peaks)


## Tables for supplement:
tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$Gene, ]
write.table(tosave.m, file="TABLES/Trendy_mouse_peakOrthologs.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$Gene, ]
write.table(tosave.h, file="TABLES/Trendy_human_peakOrthologs.csv", row.names=F, quote=F, sep=",")


## Tables for supplement:
tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[peak_genes.m[,1], ]
write.table(tosave.m, file="TABLES/Trendy_mouse_peakGenes.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[peak_genes.h[,1], ]
write.table(tosave.h, file="TABLES/Trendy_human_peakGenes.csv", row.names=F, quote=F, sep=",")



### What is the most common pattern:

# Human:
# dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up", "down", "up", "down")))
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down", "up", "down", "up"))) 
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "up", "down"))) #
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "down", "up"))) #

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "down"))) #
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "up"))) #

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down", "up"))) #
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up", "down"))) #

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up"))) #
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down"))) #




# Mouse:
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up"))) #
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down"))) #

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up", "down"))) #
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down", "up"))) #

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up", "down", "up", "down"))) #
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down", "up", "down", "up"))) #

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "up"))) #
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "down"))) #

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "up", "down"))) #




#### Step 1: Comparison of peak time in orthologs:
X <- c(subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$BreakPoint1)
Y <- c(subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$BreakPoint1)

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE, paired = TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library(ggplot2)
pdf("PLOTS/histogram_peakTime_orthologPeaks_Figure4.pdf", height=3, width=2.5)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,50), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(300), cex=.6)
dev.off()


pdf("PLOTS/spectrum_peakTime_orthologPeaks_Figure4.pdf", height=2, width=2.5)
par(mar=c(2.3,1.3,1,.1), mgp=c(1,.5,0))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', xlab="Minutes", ylab="",main="", cex.axis=.6, cex.lab=.7)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=.7, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=.7, col="cornflowerblue")
axis(2, c("Mouse", "Human"), at= c(1.7,4.2), cex=.7, tick=F)
dev.off()

# library("yarrr")
#
# X <- data.frame( Minute = X, Species = "Human")
# Y <- data.frame( Minute = Y, Species = "Mouse")
#
# longdata <- rbind(Y, X)
# 
# pdf("PLOTS/boxPlot_peakTime_orthologsPeaks_Figure4.pdf", height=7, width=5)
# par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
# pirateplot(formula = Minute ~ Species,
#            data = longdata,
#            xlab = "",inf.method = "iqr",inf.b.o = .3,point.o = .5,
#            ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
#            main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
# dev.off()


### Get peak times for all gene peaks now:

#### Step 1: Comparison of peak time in orthologs:
X <- c((peak_genes.h)$BreakPoint1)
Y <- c((peak_genes.m)$BreakPoint1)

length(X)
length(Y)

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library(ggplot2)
pdf("PLOTS/histogram_peakTime_anyPeaks_Figure4.pdf", height=3, width=2.5)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,400), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(300), cex=.6)
dev.off()



############################################################################################################################################################


# Plot specific genes

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
	       lines(tVectIn[toCol], FIT[toCol], lwd = 1.5, col=useCol)
	   }} else {
							   IDseg <- ID[1]
						       useCol <- switch(names(which.max(table(IDseg))), 
						       "0" = "black", 
						       "-1" = "cornflowerblue", 
						       "1" = "coral1")
							   	lines(tVectIn, FIT, lwd = 1.5, col=useCol)
						   }
             }

}


m1 <- extractPattern(seg.mouse, .2, Pattern = c("up", "down"), 0)
h1 <- extractPattern(seg.human, .2, Pattern = c("up", "down"), 0)
X1 <- ortho.peaks[which(ortho.peaks[,2] %in% h1$Gene & ortho.peaks[,1] %in% m1$Gene),]

## Plot first 10
toplot1 <- subset(m1, Gene %in% X1[,1])
toplot1 <- toplot1[1:8,1]
toplot1.h <- toupper(toplot1)


pdf(paste0("PLOTS/OrthologGenes_ExpressionScatter_Fig4.pdf"), height=4, width=4.7, useDingbats=F)
par(mfrow=c(2,4), mar=c(3,2,2,.5), mgp=c(1.1,.4,0))
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}
dev.off()



############################################################################################################################################################
############################################################################################################################################################







#### IER genes:
library(readxl)
ieg <- read_excel("DATA/arner_table_S5.xlsx")
ieg <- ieg[-1,4:5]
ieg <- data.frame(ieg, stringsAsFactors=F)
ieg <- ieg[1:232,]

# IEG in orthologs
ieg.h <- intersect(ieg[,1], ortho.genes.use[,2])
ieg.m <- intersect(ieg[,2], ortho.genes.use[,1])
length(ieg.h); length(ieg.m) #51

           
X = (intersect(ieg.h, ortho.peaks[,2]))
Y = (intersect(ieg.m, ortho.peaks[,1])) 
length(X)
length(Y) # 13

## IEG in -any- peaks:
ieg.h <- intersect(ieg[,1], peak_genes.h[,1]) # 34 
ieg.m <- intersect(ieg[,2], peak_genes.m[,1]) # 39 
length(ieg.h)
length(ieg.m)

# In any gene:
ieg.h <- intersect(ieg[,1], names(res.top.h$AdjustedR2)) # 78
ieg.m <- intersect(ieg[,2], names(res.top.m$AdjustedR2)) # 95
length(ieg.h)
length(ieg.m)



toplot1 <- intersect(ieg.m, ortho.peaks[,1])
toplot1.h <- intersect(ieg.h, ortho.peaks[,2])


pdf(paste0("PLOTS/OrthologGenes_ExpressionScatter_IEG_ALL.pdf"), height=4, width=8, useDingbats=F)
par(mfrow=c(2,4), mar=c(3,2,2,.5), mgp=c(1.1,.4,0))
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}

for(i in 5:8){
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 5:8){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}


for(i in 9:12){
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 9:12){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}

for(i in 13){3
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 13){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}
dev.off()







toplot1 <- c("Myc", "Sqstm1", "Adm", "Ubc")
toplot1.h <- toupper(toplot1)


pdf(paste0("PLOTS/OrthologGenes_ExpressionScatter_IEG_Fig4.pdf"), height=4, width=4.7, useDingbats=F)
par(mfrow=c(2,4), mar=c(3,2,2,.5), mgp=c(1.1,.4,0))
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.m, tVectIn=t.v.m,
            featureNames = toplot1[i], 
            trendyOutData = seg.mouse)
}
for(i in 1:4){
XX<- fancyPlot(data.norm.scale.h, tVectIn=t.v.h,
            featureNames = toplot1.h[i], 
            trendyOutData = seg.human)
}
dev.off()

