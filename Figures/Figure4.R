setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)

### Get Orthologs:

# Clean this up a bit for future use:
orth.genes.clean1 <- orth.genes[which(orth.genes[,1]!="" & orth.genes[,2]!=""),]
orth.genes.clean <- orth.genes.clean1[!duplicated(orth.genes.clean1),]

# For using a cutoff of .5:
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

top1 <- merge(orth.genes.clean, top.human, by="hgnc_symbol")
top2 <- merge(top1,top.mouse, by="mgi_symbol")
top2 <- top2[!duplicated(top2),]
dupg <- top2[which(duplicated(top2[,2])),2] 
subset(top2, hgnc_symbol %in% dupg)
dupg <- top2[which(duplicated(top2[,1])),1] 
subset(top2, mgi_symbol %in% dupg)
TORM <- c(45, 32)
top2 <- top2[-TORM,]
dim(top2)

# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse, !(Gene %in% top2$mgi_symbol))
intersect(toupper(mouse.genes.check1[,1]), top.human$Gene)
recovered.g <- data.frame(mgi_symbol = c("Ubc", "Iffo2"), hgnc_symbol = c("UBC","IFFO2"), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)

ortho.genes.use <- top2

# Ortholog Peaks:
								
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


### What is the most common pattern:

# Human:
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up", "down", "up", "down"))) #0
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down", "up", "down", "up"))) #0
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "up", "down"))) #6
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "down", "up"))) #0

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "down"))) #3
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("same", "up"))) #9

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down", "up"))) #14
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up", "down"))) #33

dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("up"))) #5
dim(extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down"))) #24




# Mouse:
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up"))) #3
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down"))) #2

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up", "down"))) #57
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down", "up"))) #24

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("up", "down", "up", "down"))) #5
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down", "up", "down", "up"))) #1

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "up"))) #13
dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "down"))) #6

dim(extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("same", "up", "down"))) #7




#### Step 1: Comparison of peak time in orthologs:
X <- c(subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)$BreakPoint1)
Y <- c(subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)$BreakPoint1)


PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
PP
library(ggplot2)
pdf("PLOTS/histogram_peakTime_orthologPeaks_Figure4.pdf", height=8, width=12)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,15), border="brown3", 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=2.5, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6),border="dodgerblue3",  breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6),alpha("brown1", .6)), cex=2)
dev.off()

pdf("PLOTS/spectrum_peakTime_orthologPeaks_Figure4.pdf", height=3.5, width=16)
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

pdf("PLOTS/boxPlot_peakTime_orthologsPeaks_Figure4.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
pirateplot(formula = Minute ~ Species,
           data = longdata,
           xlab = "",inf.method = "iqr",inf.b.o = .3,point.o = .5,
           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()


# Average time: 
peak_genes.h.ortho <- subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)
rownames(peak_genes.h.ortho) <- peak_genes.h.ortho$Gene
peak_genes.m.ortho <- subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)
rownames(peak_genes.m.ortho) <- peak_genes.m.ortho$Gene

X <- as.vector(peak_genes.h.ortho$BreakPoint1)
names(X) <- rownames(peak_genes.h.ortho)

Y <- as.vector(peak_genes.m.ortho$BreakPoint1)
names(Y) <- rownames(peak_genes.m.ortho)


mean(Y[ortho.peaks$mgi_symbol] - X[ortho.peaks$hgnc_symbol], na.rm=T)

### 278 minutes earlier


## Not enough genes with dip or immediate down ortholog to analyze them as a group.



### Get peak times for all gene peaks now:

#### Step 1: Comparison of peak time in orthologs:
X <- c(subset(peak_genes.h)$BreakPoint1)
Y <- c(subset(peak_genes.m)$BreakPoint1)


PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
PP
library(ggplot2)
pdf("PLOTS/histogram_peakTime_anyPeaks_Figure4.pdf", height=8, width=12)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,80), border="brown3", 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=2.5, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6),border="dodgerblue3",  breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6),alpha("brown1", .6)), cex=2)
dev.off()

pdf("PLOTS/spectrum_peakTime_anyPeaks_Figure4.pdf", height=3.5, width=16)
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

pdf("PLOTS/boxPlot_peakTime_anyPeaks_Figure4.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
pirateplot(formula = Minute ~ Species,
           data = longdata,
           xlab = "",inf.method = "iqr",inf.b.o = .3,point.o = .5,
           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()


#
#
# ## Tried but not enough data to do more analysis. Focus on Immediate Down and Up later.
# #### Step 1: Comparison of dip time in orthologs:
#
# dip_genes.m <- extractPattern(seg.mouse[ortho.genes.use[,1]], Pattern =c("down", "same"))
# dip_genes.h <- extractPattern(seg.human[ortho.genes.use[,2]], Pattern =c("down", "same"))
#
# dip_genes.m <- dip_genes.m[which(!duplicated(dip_genes.m[,1])),]
# dip_genes.h <- dip_genes.h[which(!duplicated(dip_genes.h[,1])),]
#
# set1 <- subset(ortho.genes.use, hgnc_symbol %in% dip_genes.h$Gene)
# ortho.dips <- subset(set1, mgi_symbol %in% dip_genes.m$Gene)
#
# X <- c(subset(dip_genes.h, Gene %in% ortho.dips$hgnc_symbol)$BreakPoint1)
# Y <- c(subset(dip_genes.m, Gene %in% ortho.dips$mgi_symbol)$BreakPoint1)
#
#
# PP <- round(wilcox.test(X,Y)$p.value, 3)
# if( PP < .001) {PP <- "< .001"}
# PP
# library(ggplot2)
# par(mar=c(6,6,3,1), mgp=c(4,1,0))
# hist(X, xlim=c(0,600), ylim=c(0,15), border="brown3",
#   col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20),
#   main="", xlab="Minute",
#   cex.axis=2.5, cex.lab=3)
# hist(Y, add=T,  col=alpha("cornflowerblue", .6),border="dodgerblue3",  breaks = seq(0, 600, length.out=20))
# legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6),alpha("brown1", .6)), cex=2)
# # dev.off()
#
# par(mar=c(3,3,1,1), mfrow=c(1,1), mgp=c(0,2,0))
# plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', ylab="", xlab="",main="", cex.axis=3)
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
# points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=2, col="brown1", bg="black")
# points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=2, col="cornflowerblue")
# # dev.off()
#
# library("yarrr")
#
# X <- data.frame( Minute = X, Species = "Human")
# Y <- data.frame( Minute = Y, Species = "Mouse")
#
# longdata <- rbind(Y, X)
#
# par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
# pirateplot(formula = Minute ~ Species,
#            data = longdata,
#            xlab = "",inf.method = "iqr",inf.b.o = .3,point.o = .5,
#            ylab = "Minutes", pal=c("cornflowerblue", "brown1"),
#            main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
# # dev.off()





#### IER genes:

ieg <- read.csv("~/Downloads/arner_table_S5.csv", header=T, stringsAsFactors=F)
ieg <- ieg[-1,]
ieg <- ieg[1:233,4:5]


# IEG in orthologs

ieg.h <- intersect(ieg[,1], ortho.genes.use[,2])
ieg.m <- intersect(ieg[,2], ortho.genes.use[,1])
length(ieg.h); length(ieg.m) #1 0


# apply(res.top.h$Segment.Trends[ieg.h,], 1, function(x) {
#         x[x == -1] <- 2
#         if (length(x) == 1) {
#             x = rep(x, 1)
#         }
#         x <- paste(na.omit(x), collapse = "")
#         return(x)
#     })



intersect(ieg.h, ortho.peaks[,2])
intersect(ieg.m, ortho.peaks[,1]) # 4





## IEG in -any- peaks:

ieg.h <- intersect(ieg[,1], peak_genes.h[,1]) # 16 / 121
ieg.m <- intersect(ieg[,2], peak_genes.m[,1]) #15 / 362
length(ieg.h)
length(ieg.m)

# In any gene:
ieg.h <- intersect(ieg[,1], names(res.top.h$AdjustedR2)) # 23 / 739
ieg.m <- intersect(ieg[,2], names(res.top.m$AdjustedR2)) #37 / 749
length(ieg.h)
length(ieg.m)


sum(res.top.m$Segment.Trends[ieg.m,1] == -1)
sum(res.top.m$Segment.Trends[ieg.m,1] == 1)
sum(res.top.m$Segment.Trends[ieg.m,1] == 0)


sum(res.top.h$Segment.Trends[ieg.h,1] == -1)
sum(res.top.h$Segment.Trends[ieg.h,1] == 1)
sum(res.top.h$Segment.Trends[ieg.h,1] == 0)


#
# ## Not interesting, skip for now.
#
# use.genes1 <- rbind(extractPattern(seg.mouse, Pattern =c("down", "up")), extractPattern(seg.mouse, Pattern =c("down", "same")))
# use.genes2 <- extractPattern(seg.mouse, Pattern =c("up", "down"))
# use.genes <- unique(setdiff(use.genes1$Gene, use.genes2$Genes))
# ieg.m <- intersect(ieg[,2], use.genes) #
# length(ieg.m)
#
#
# use.genes1 <- rbind(extractPattern(seg.human, Pattern =c("down", "up")), extractPattern(seg.human, Pattern =c("down", "same")))
# use.genes1 <- rbind(extractPattern(seg.human, Pattern =c("same", "down")), extractPattern(seg.human, Pattern =c("down", "same")))
# use.genes2 <- extractPattern(seg.human, Pattern =c("up", "down"))
# use.genes <- unique(setdiff(use.genes1$Gene, use.genes2$Genes))
# ieg.h <- intersect(ieg[,1], use.genes) #
# length(ieg.h)
#
#
# use.g <- intersect(ieg[,1], rownames(res.top.h$Segment.Trends) )
#
# apply(res.top.h$Segment.Trends[use.g,], 1, function(x) {
#         x[x == -1] <- 2
#         if (length(x) == 1) {
#             x = rep(x, 1)
#         }
#         x <- paste(na.omit(x), collapse = "")
#         return(x)
#     })




############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################




# Plot specific genes

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


m1 <- extractPattern(seg.mouse, .5, Pattern = c("up", "down"), 0)
h1 <- extractPattern(seg.human, .5, Pattern = c("up", "down"), 0)
X1 <- ortho.peaks[which(ortho.peaks[,2] %in% h1$Gene & ortho.peaks[,1] %in% m1$Gene),]

## Plot first 10
toplot1 <- subset(m1, Gene %in% X1[,1])
toplot1 <- toplot1[1:10,1]
toplot1.h <- toupper(toplot1)


toplot1 <- c("Myc", "Sqstm1", "Adm", "Ubc")
toplot1.h <- toupper(toplot1)

for(i in 1:4){

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
