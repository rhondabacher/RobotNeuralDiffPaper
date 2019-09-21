setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

# Peak genes, divide into initial slope direction:
library(Trendy)
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)
dim(ortho.peaks)

## For peak genes get the UP slope and DOWN slope
peak.com.h <- subset(peak_genes.h, Gene %in% ortho.peaks$hgnc_symbol)
whichTime.up.h<-c()
for(i in 1:nrow(peak.com.h)) {
	keep <- which(round(res.top.h$Breakpoints[peak.com.h[i,1],]) == round(c(peak.com.h[i,2])))
	if(keep == 1) {
		whichTime.up.h[i] <- 0
	} else {
	whichTime.up.h[i] <- res.top.h$Breakpoints[peak.com.h[i,1],(keep-1)]
	}
	
}
names(whichTime.up.h) <- peak.com.h[,1]

peak.com.m <- subset(peak_genes.m, Gene %in% ortho.peaks$mgi_symbol)
whichTime.up.m<-c()
for(i in 1:nrow(peak.com.m)) {
	keep <- which(round(res.top.m$Breakpoints[peak.com.m[i,1],]) == round(c(peak.com.m[i,2])))
	if(keep == 1) {
		whichTime.up.m[i] <- 0
	} else {
	whichTime.up.m[i] <- res.top.m$Breakpoints[peak.com.m[i,1],(keep-1)]
	}
	
}
names(whichTime.up.m) <- peak.com.m[,1]  

pcntStart0.h <- mean(whichTime.up.h==0)*100
pcntStart0.m <- mean(whichTime.up.m==0)*100


propTestVals <- prop.test(c(sum(whichTime.up.h==0), sum(whichTime.up.m==0)), c(length(whichTime.up.h), length(whichTime.up.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/percent_FirstTime_Up_CommonPeaks_Figure5.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(sum(whichTime.up.m == 0), sum(whichTime.up.h==0), sum(whichTime.up.m > 0), sum(whichTime.up.h > 0)),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,150), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()



############################################################################################
## For ANY peak genes get the UP slope and DOWN slope
peak.com.h <- peak_genes.h
whichTime.up.h<-c()
for(i in 1:nrow(peak.com.h)) {
	keep <- which(round(res.top.h$Breakpoints[peak.com.h[i,1],]) == round(c(peak.com.h[i,2])))
	if(keep == 1) {
		whichTime.up.h[i] <- 0
	} else {
	whichTime.up.h[i] <- res.top.h$Breakpoints[peak.com.h[i,1],(keep-1)]
	}
	
}
names(whichTime.up.h) <- peak.com.h[,1] 


peak.com.m <- peak_genes.m
whichTime.up.m<-c()
for(i in 1:nrow(peak.com.m)) {
	keep <- which(round(res.top.m$Breakpoints[peak.com.m[i,1],]) == round(c(peak.com.m[i,2])))
	if(keep == 1) {
		whichTime.up.m[i] <- 0
	} else {
	whichTime.up.m[i] <- res.top.m$Breakpoints[peak.com.m[i,1],(keep-1)]
	}
	
}
names(whichTime.up.m) <- peak.com.m[,1]

pcntStart0.h <- mean(whichTime.up.h==0)*100
pcntStart0.m <- mean(whichTime.up.m==0)*100

# Obviously
propTestVals <- prop.test(c(sum(whichTime.up.h==0), sum(whichTime.up.m==0)), c(length(whichTime.up.h), length(whichTime.up.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/percent_FirstTime_Up_AnyPeaks_Figure5.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(pcntStart0.m, pcntStart0.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="% Genes", xlab="",
ylim=c(0,100), cex.axis=.6, cex.lab=.7
)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()






############################################################################################

## Get common UP

# Above was peak orthologs, now just genes going UP and orthologs
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

upgenes.m <- c(names(which(res.top.m$Segment.Trends[,1] == 0 & res.top.m$Segment.Trends[,2] == 1)), 
                names(which(res.top.m$Segment.Trends[,1] == 1)))

upgenes.h <- c(names(which(res.top.h$Segment.Trends[,1] == 0 & res.top.h$Segment.Trends[,2] == 1)), 
                names(which(res.top.h$Segment.Trends[,1] == 1)))

## Genes in common between mouse and human # Orthologs

## Orthologs going Up:
set1 <- subset(ortho.genes.use, hgnc_symbol %in% upgenes.h)
ortho.up <- subset(set1, mgi_symbol %in% upgenes.m)


# Get time of UP or delayed UP
notzero <- names(which(res.top.m$Segment.Trends[ortho.up[,1],1] == 0 & res.top.m$Segment.Trends[ortho.up[,1],2] == 1))
upbp.m <- res.top.m$Breakpoints[notzero,1]
upbp.m <- c(upbp.m, rep(0, length(which(res.top.m$Segment.Trends[ortho.up[,1],1] == 1))))

notzero <- names(which(res.top.h$Segment.Trends[ortho.up[,2],1] == 0 & res.top.h$Segment.Trends[ortho.up[,2],2] == 1))
upbp.h <- res.top.h$Breakpoints[notzero,1]
upbp.h <- c(upbp.h, rep(0, length(which(res.top.h$Segment.Trends[ortho.up[,2],1] == 1))))

pcntStart0.h <- mean(upbp.h==0)*100
pcntStart0.m <- mean(upbp.m==0)*100


propTestVals <- prop.test(c(sum(upbp.h==0), sum(upbp.m==0)), c(length(upbp.h), length(upbp.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/percent_FirstTime_Up_CommonUP_Figure5.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(sum(upbp.m == 0), sum(upbp.h==0), sum(upbp.m > 0), sum(upbp.h > 0)),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,300), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()



############### ALL UP GENES:
# Above was up orthologs, now just genes going UP


library(Trendy)
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

upgenes.m <- c(names(which(res.top.m$Segment.Trends[,1] == 0 & res.top.m$Segment.Trends[,2] == 1)), 
                names(which(res.top.m$Segment.Trends[,1] == 1)))

upgenes.h <- c(names(which(res.top.h$Segment.Trends[,1] == 0 & res.top.h$Segment.Trends[,2] == 1)), 
                names(which(res.top.h$Segment.Trends[,1] == 1)))

notzero <- names(which(res.top.m$Segment.Trends[upgenes.m,1] == 0 & res.top.m$Segment.Trends[upgenes.m,2] == 1))
upbp.m <- res.top.m$Breakpoints[notzero,1]
upbp.m <- c(upbp.m, rep(0, length(which(res.top.m$Segment.Trends[upgenes.m,1] == 1))))

notzero <- names(which(res.top.h$Segment.Trends[upgenes.h,1] == 0 & res.top.h$Segment.Trends[upgenes.h,2] == 1))
upbp.h <- res.top.h$Breakpoints[notzero,1]
upbp.h <- c(upbp.h, rep(0, length(which(res.top.h$Segment.Trends[upgenes.h,1] == 1))))


pcntStart0.h <- mean(upbp.h==0)*100
pcntStart0.m <- mean(upbp.m==0)*100

# Obviously
propTestVals <- prop.test(c(sum(upbp.h==0), sum(upbp.m==0)), c(length(upbp.h), length(upbp.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/percent_FirstTime_Up_allUP_Figure5.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(pcntStart0.m, pcntStart0.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="% Genes", xlab="",
ylim=c(0,100), cex.axis=.6, cex.lab=.7
)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
# legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()



