setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)


# Dynamic genes:
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)


# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends, 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends, 1, function(x) which(x != 0)[1])]

sum(res.top.m$Trends[,1] !=0) / nrow(res.top.m$Trends) # 91%
sum(res.top.m$Trends[,1] !=0) / nrow(res.top.m$Trends) # 93%

X <- timeUpDown.h
Y <- timeUpDown.m

library(ggplot2)

pdf("PLOTS/histogram_FirstTime_UpDown_Figure2.pdf", height=6, width=8)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,700), 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=50), 
	main="", xlab="Minute",
	cex.axis=2, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), breaks = seq(0, 600, length.out=50))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2, bty='n')
dev.off()


## Now get what percent of genes go up/down at time 0 versus delayed up/down
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[,], 1, function(x) which(x != 0)[1])]

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100

prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)))
prop.test(c(sum(timeUpDown.h!=0), sum(timeUpDown.m!=0)), c(length(timeUpDown.h), length(timeUpDown.m)))

pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_Figure2.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()


## Monotonic genes, only up or down the entire time:

onlyDown.m <- sum(apply(res.top.m$Trends, 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends, 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends, 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends, 1, function(x) all(x == 1)))

prop.test(c((onlyDown.m), (onlyDown.h)), c(length(res.top.m$Trends), length(res.top.h$Trends)))
prop.test(c((onlyUp.m), (onlyUp.h)), c(length(res.top.m$Trends), length(res.top.h$Trends)))


pdf("PLOTS/numGenes_Monotonic_UpDown_allTrendyGenes_Figure2.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(onlyUp.m, onlyUp.h, onlyDown.m, onlyDown.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,250), cex.axis=2
)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()


### Now histogram of the END of the up/down segments. When does that initial trend change?

pdf("PLOTS/histogram_trendUntilTime_UpDown_Figure2.pdf", height=7, width=18)
par(mar=c(5.5,5,1,1), mfrow=c(1,2), mgp=c(3.5,1,0))

res.top.m$Breakpoints[,1][is.na(res.top.m$Breakpoints[,1])] <- 600
res.top.h$Breakpoints[,1][is.na(res.top.h$Breakpoints[,1])] <- 600

#Up
up.bp.m <- round(c(res.top.m$Breakpoints[names(which(res.top.m$Segment.Trends[,1] == 1)),1]))
up.bp.h <- round(c(res.top.h$Breakpoints[names(which(res.top.h$Segment.Trends[,1] == 1)),1]))


hist(up.bp.h, xlim=c(0,600), ylim=c(0,140), 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=100), border="brown3", 
	main="", xlab="Minute", ylab="",
	cex.axis=2.5, cex.lab=3)
hist(up.bp.m, 
	col=alpha("cornflowerblue", .6), breaks = seq(0, 600, length.out=100), border="dodgerblue3", 
	add=TRUE)
legend('topleft', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2, bty='n')


## Down
down.bp.m <- round(c(res.top.m$Breakpoints[names(which(res.top.m$Segment.Trends[,1] == -1)),1]))
down.bp.h <- round(c(res.top.h$Breakpoints[names(which(res.top.h$Segment.Trends[,1] == -1)),1]))


hist(down.bp.h, xlim=c(0,600), ylim=c(0,201), 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=100), border="brown3", 
	main="", xlab="Minute", ylab="",
	cex.axis=2.5, cex.lab=3)
hist(down.bp.m, 
	col=alpha("cornflowerblue", .6), breaks = seq(0, 600, length.out=100), border="dodgerblue3", 
	add=TRUE)
legend('topleft', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2, bty='n')

dev.off()

