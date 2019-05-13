setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")


# Peak genes, divide into initial slope direction:
library(Trendy)
peak_genes.m <- peak_genes.m[which(!duplicated(peak_genes.m[,1])),]
peak_genes.h <- peak_genes.h[which(!duplicated(peak_genes.h[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m$Gene)


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

# Obviously
prop.test(c(sum(whichTime.up.h==0), sum(whichTime.up.m==0)), c(length(whichTime.up.h), length(whichTime.up.m)))



pdf("PLOTS/percent_FirstTime_Up_CommonPeaks_Figure5.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()


# Save these for later plot:

whichTime.up.h.orthopeaks <- whichTime.up.h
whichTime.up.m.orthopeaks <- whichTime.up.m








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
prop.test(c(sum(whichTime.up.h==0), sum(whichTime.up.m==0)), c(length(whichTime.up.h), length(whichTime.up.m)))



pdf("PLOTS/percent_FirstTime_Up_AnyPeaks_Figure5.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()


# Save these for later plot:

whichTime.up.h.anypeaks <- whichTime.up.h
whichTime.up.m.anypeaks <- whichTime.up.m













## Get common UP

# Above was peak orthologs, now just genes going UP and orthologs
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

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


upTime_commonUP.h <- upbp.h
upTime_commonUP.m <- upbp.m

prop.test(c(sum(upTime_commonUP.h==0), sum(upTime_commonUP.m==0)), 
                    c(length(upTime_commonUP.h), length(upTime_commonUP.m)))

pdf("PLOTS/percent_FirstTime_Up_CommonUP_Figure5.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()


############### ALL UP GENES:
# Above was up orthologs, now just genes going UP
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

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

prop.test(c(sum(upbp.h==0), sum(upbp.m==0)), 
                    c(length(upbp.h), length(upbp.m)))

pdf("PLOTS/percent_FirstTime_Up_allUP_Figure5.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()

upTime_allUP.h <- upbp.h
upTime_allUP.m <- upbp.m



### Plot together

pdf("PLOTS/percent_MeanOf_FirstTime_Up_Figure5.pdf", height=5, width=18)
par(mar=c(3,3,3,1))
barplot(c(mean(upTime_allUP.m), mean(upTime_allUP.h),
					mean(upTime_commonUP.m), mean(upTime_commonUP.h), 
					mean(whichTime.up.m.orthopeaks), mean(whichTime.up.h.orthopeaks),
					mean(whichTime.up.m.anypeaks), mean(whichTime.up.h.anypeaks))), 
space=c(.5,.1,1,.1, 1, .1, 1, .1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,600), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()




