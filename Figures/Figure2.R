setwd("~/RobotSeq/")



load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)


res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)


## Tables for supplement:

# All dynamic genes in each species
tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[names(sort(res.top.m$Breakpoints[,1],na.last=TRUE)), ]
write.table(tosave.m, file="TABLES/Trendy_mouse_AllDynamic.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[names(sort(res.top.h$Breakpoints[,1],na.last=TRUE)), ]
write.table(tosave.h, file="TABLES/Trendy_mouse_AllDynamic.csv", row.names=F, quote=F, sep=",")

dim(tosave.m) # 3721
dim(tosave.h) # 4332

# Gene dynamic in ONLY one species
tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[names(sort(res.top.m$Breakpoints[setdiff(rownames(tosave.m), ortho.genes.use$mgi_symbol),1],na.last=TRUE)), ]
dim(tosave.m) # 2285
write.table(tosave.m, file="TABLES/Trendy_mouseOnly_Dynamic.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[names(sort(res.top.h$Breakpoints[setdiff(rownames(tosave.h), ortho.genes.use$hgnc_symbol),1],na.last=TRUE)), ]
dim(tosave.h) # 2896
write.table(tosave.h, file="TABLES/Trendy_humanOnly_Dynamic.csv", row.names=F, quote=F, sep=",")


tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[names(sort(res.top.m$Breakpoints[ortho.genes.use$mgi_symbol,1],na.last=TRUE)), ]
dim(tosave.m) # 1436
write.table(tosave.m, file="TABLES/Trendy_mouse_OrthologDynamic.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[names(sort(res.top.h$Breakpoints[ortho.genes.use$hgnc_symbol,1],na.last=TRUE)), ]
dim(tosave.h) # 1436
write.table(tosave.h, file="TABLES/Trendy_human_OrthologDynamic.csv", row.names=F, quote=F, sep=",")


#### Now analyze orthologs

# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) which(x != 0)[1])]

library(ggplot2)

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100

propTestVals <-prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 2)

pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_Figure2_Orthologs.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(sum(timeUpDown.h == 0), sum(timeUpDown.m==0), sum(timeUpDown.h > 0), sum(timeUpDown.m > 0)),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,1500), cex.axis=.6, cex.lab=.7
)
mtext(c("Minute = 0", "Minute > 0"), side=1, at = c(1.5, 4.5), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(2.5), cex=.6)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
dev.off()




## ORTHOLOGS Monotonic genes, only up or down the entire time:

onlyDown.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == 1)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

propTestVals <- prop.test(c((all.mono.h), (all.mono.m)), c(length(ortho.genes.use$mgi_symbol), length(ortho.genes.use$mgi_symbol)), conf.level=.99)
propTestVals <- round(propTestVals$conf.int[1:2]*100, 2)

pdf("PLOTS/numGenes_Monotonic_allTrendyGenes_Figure2_Orthologs.pdf", height=2, width=2)
par(mar=c(1.5,3,1,.1), mgp=c(2,1,0))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="", ylab="# Orthologs", xlab="",
ylim=c(0,850), cex.axis=.6, cex.lab=.7
)
mtext(c("Mouse", "Human"), side=1, at = c(.5,1.8), cex=.6)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=.5, bty='n')
mtext(bquote("99% CI " ~ Delta~ "P%: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1), cex=.6)
dev.off()

all.mono.h
all.mono.m
all.mono.h / all.mono.m


########################################################################################################
## Get breakpoints and slopes for ortholog genes

library(gplots)
library(ggplot2)

all.slopes.mouse <- res.top.m$Segment.Slopes
colnames(all.slopes.mouse) <- paste0("mouse_slope", seq_len(ncol(all.slopes.mouse)))

all.bp.mouse <- res.top.m$Breakpoints
colnames(all.bp.mouse) <- paste0("mouse_breakpoint", seq_len(ncol(all.bp.mouse)))

all.slopes.human <- res.top.h$Segment.Slopes
colnames(all.slopes.human) <- paste0("human_slope", seq_len(ncol(all.slopes.human)))

all.bp.human <- res.top.h$Breakpoints
colnames(all.bp.human) <- paste0("human_breakpoint", seq_len(ncol(all.bp.human)))


all.bp.human <- (all.bp.human[ortho.genes.use$hgnc_symbol,])
all.slopes.human <- (all.slopes.human[ortho.genes.use$hgnc_symbol,])
all.bp.mouse <- (all.bp.mouse[ortho.genes.use$mgi_symbol,])
all.slopes.mouse <- (all.slopes.mouse[ortho.genes.use$mgi_symbol,])

X <- (c(all.bp.human[,1]))
Y <- (c(all.bp.mouse[,1]))

sum(X <= 250, na.rm=T) / length(X[!is.na(X)])
sum(Y <= 250, na.rm=T) / length(Y[!is.na(Y)])

propTestVals <- wilcox.test(X, Y, conf.level = .99, conf.int=TRUE, paired = TRUE)
propTestVals <- round(propTestVals$conf.int[1:2], 0)


pdf("PLOTS/histogram_firstBreakpoint_Orthologs_Figure2.pdf", height=2, width=3)
par(mar=c(2.1,2,1,.1), mgp=c(1.2,.5,0))
hist(X, xlim=c(0,600), ylim=c(0,300), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=.6, cex.lab=.7)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=.5)
dev.off()

pdf("PLOTS/spectrum_firstBreakpoint_Orthologs_Figure2.pdf", height=2, width=4.2)
par(mar=c(2.3,1.5,1,.1), mgp=c(1,.5,0))
plot(X, rep(2, length(X)), ylim=c(1,5), xlim=c(0,600), yaxt='n', xlab="Minutes", ylab="",main="", cex.axis=.6, cex.lab=.7)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
points(X, rep(2, length(X)), ylim=c(0,10), xlim=c(0,600), pch=250, font=5, cex=.7, col="brown1", bg="black")
points(Y, rep(4, length(Y)), ylim=c(0,10), pch=250, font=5, cex=.7, col="cornflowerblue")
axis(2, c("Mouse", "Human"), at= c(1.7,4.2), cex=.7, tick=F)
dev.off()

library("yarrr")

X <- data.frame( Minute = X, Species = "Human")
Y <- data.frame( Minute = Y, Species = "Mouse")

longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_firstBreakpoint_Orthologs_Figure2.pdf", height=2, width=2.5)
par(mar=c(1,2.5,1,.1), mgp=c(1.6,.5,0))
pirateplot(formula = Minute ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "M: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
