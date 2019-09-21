setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

library(Trendy)

###################################################################################################
###################################################################################################
###################################################################################################
## Not using orthologs here!

## Breaks before various times:

res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

sum(res.top.m$Breakpoint <= 60, na.rm=T) # 504
length(unique(c(names(which(res.top.m$Breakpoint[,3] <= 60)),
        names(which(res.top.m$Breakpoint[,2] <= 60)), names(which(res.top.m$Breakpoint[,1] <= 60))))) # 497
sum(res.top.m$Breakpoint <= 100, na.rm=T) # 892
sum(res.top.m$Breakpoint <= 250, na.rm=T) # 2919

sum(res.top.h$Breakpoint <= 60, na.rm=T) # 29
length(unique(c(names(which(res.top.h$Breakpoint[,3] <= 60)),
        names(which(res.top.h$Breakpoint[,2] <= 60)), names(which(res.top.h$Breakpoint[,1] <= 60))))) # 29
sum(res.top.h$Breakpoint <= 100, na.rm=T) # 138
sum(res.top.h$Breakpoint <= 250, na.rm=T) # 457

X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)

library("yarrr")

X <- data.frame( Breakpoint = X, Species = "Human")
Y <- data.frame( Breakpoint = Y, Species = "Mouse")
longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_AllBreakpoints_anyGenes.pdf", height=3, width=2)
par(mar=c(1,2.5,1,.1), mgp=c(1.6,.5,0))
pirateplot(formula = Breakpoint ~ Species,
	           data = longdata,
	           xlab = "", inf.b.o = .3,point.o = .5, 
	           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3, bar.lwd=1, cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()


X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))


pdf("PLOTS/histogram_GenesInFig1_AllBreakpoints.pdf", height=4, width=6)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0), mfrow=c(2,1))
hist(Y, xlim=c(0,600), ylim=c(0,500), 
	col="cornflowerblue", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="Mouse", xlab="Minute",border="dodgerblue3", 
	cex.axis=.6, cex.lab=.7, cex.main=.8)
axis(1, at=c(60, 250), cex.axis=.6, cex.lab=.7)
abline(v = c(60, 100, 250), lwd=1, lty=2, col="gray60")
hist(X, xlim=c(0,600), ylim=c(0,400), main="Human", cex.main=.8,
	col="brown1", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", xlab="Minute",border="brown3", 
	cex.axis=.6, cex.lab=.7)
axis(1, at=c(60, 250), cex.axis=.6, cex.lab=.7)
abline(v = c(60, 100, 250), lwd=1, lty=2, col="gray60")
dev.off()







# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.m <- timeUpDown.m[!is.na(timeUpDown.m)]
timeUpDown.h <- timeUpDown.h[!is.na(timeUpDown.h)]
sum(res.top.h$Trends[,1] !=0) / length(timeUpDown.h) # 92%
sum(res.top.m$Trends[,1] !=0) / length(timeUpDown.m) # 94%

propTestVals <- prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100


pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_Figure3.pdf", height=3, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(mean(timeUpDown.h == 0), mean(timeUpDown.m==0), mean(timeUpDown.h > 0), mean(timeUpDown.m > 0))*100,
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="% Genes", xlab="",
ylim=c(0,100), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()


 
onlyDown.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == 1)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

propTestVals <- prop.test(c((all.mono.h), (all.mono.m)), c(nrow(res.top.h$Trends), nrow(res.top.m$Trends)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 3)

pdf("PLOTS/numGenes_Monotonic_allTrendyGenes_Figure.pdf", height=3, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,2550), cex.axis=.6, cex.lab=.7
)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()

all.mono.h
all.mono.m
all.mono.h / all.mono.m