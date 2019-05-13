setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)


res.top <- topTrendy(seg.mouse, .5)

# How many breakpoints happen before 60 minutes?
sum(res.top$Breakpoint <= 60, na.rm=T)
length(unique(c(names(which(res.top$Breakpoint[,2] <= 60)), names(which(res.top$Breakpoint[,1] <= 60)))))

sum(res.top$Breakpoint <= 100, na.rm=T) #264
sum(res.top$Breakpoint <= 250, na.rm=T) #738


# Human

# All genes with an adjusted R^2 >= .5
res.top <- topTrendy(seg.human, .5)

# How many breakpoints happen before 30 minutes?
sum(res.top$Breakpoint[,1] <= 60, na.rm=T)
# 0
res.top$Breakpoint[names(which(res.top$Breakpoint[,2] <= 60)),]

length(unique(c(names(which(res.top$Breakpoint[,2] <= 60)), names(which(res.top$Breakpoint[,1] <= 60)))))


sum(res.top$Breakpoint <= 100, na.rm=T) # 38
sum(res.top$Breakpoint <= 250, na.rm=T) # 100




###################################################################################################
###################################################################################################
###################################################################################################

## Joint plots

# Dynamic genes in common:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

sum(res.top.m$Breakpoint <= 60, na.rm=T)
sum(res.top.m$Breakpoint <= 100, na.rm=T)
sum(res.top.m$Breakpoint <= 250, na.rm=T)

sum(res.top.h$Breakpoint <= 60, na.rm=T)
sum(res.top.h$Breakpoint <= 100, na.rm=T)
sum(res.top.h$Breakpoint <= 250, na.rm=T)

X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))

wilcox.test(X,Y)
PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
PP

library("yarrr")

X <- data.frame( Breakpoint = X, Species = "Human")
Y <- data.frame( Breakpoint = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_GenesInFig1_AllBreakpoints.pdf", height=12, width=8)
par(mar=c(4,6,2,1), mgp = c(4, 1, 0))
pirateplot(formula = Breakpoint ~ Species,
           data = longdata,
           xlab = "", inf.method = "iqr",inf.b.o = .3,point.o = .5,
           ylab = "Minute", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=3, cex.axis=2,cex.names=3)
dev.off()

X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))


pdf("PLOTS/histogram_GenesInFig1_AllBreakpoints.pdf", height=12, width=15)
par(mar=c(5,5,4,2), mfrow=c(2,1))
hist(Y, xlim=c(0,600), ylim=c(0,80), 
	col="cornflowerblue", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="", xlab="",border="dodgerblue3", 
	cex.axis=2, cex.lab=3)
axis(1, at=c(60, 250), cex.axis=2, cex.lab=3)
abline(v = c(60, 100, 250), lwd=1, lty=2, col="gray60")
hist(X, xlim=c(0,600), ylim=c(0,80), 
	col="brown1", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="", xlab="",border="brown3", 
	cex.axis=2, cex.lab=3)
axis(1, at=c(60, 250), cex.axis=2, cex.lab=3)
abline(v = c(60, 100, 250), lwd=1, lty=2, col="gray60")
dev.off()







# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[,], 1, function(x) which(x != 0)[1])]

sum(res.top.h$Trends[,1] !=0) / length(timeUpDown.h) # 93%
sum(res.top.m$Trends[,1] !=0) / length(timeUpDown.m) # 91%

prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)))
prop.test(c(sum(timeUpDown.h!=0), sum(timeUpDown.m!=0)), c(length(timeUpDown.h), length(timeUpDown.m)))


pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100


pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_Figure2.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()






onlyDown.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[,], 1, function(x) all(x == 1)))

prop.test(c((onlyDown.m), (onlyDown.h)), c(nrow(res.top.m$Trends), nrow(res.top.h$Trends)))
prop.test(c((onlyUp.m), (onlyUp.h)), c(nrow(res.top.m$Trends), nrow(res.top.h$Trends)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

prop.test(c((all.mono.m), (all.mono.h)), c(nrow(res.top.m$Trends), nrow(res.top.h$Trends)))

pdf("PLOTS/numGenes_Monotonic_allTrendyGenes_Figure2.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,550), cex.axis=2
)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()

