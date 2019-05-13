setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)



# For Ortholog Genes, what % of genes have an iniital trend?
# Dynamic genes:
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

## Tables for supplement:

tosave.m <- formatResults(res.top.m)
tosave.m <- tosave.m[names(sort(res.top.m$Breakpoints[ortho.genes.use$mgi_symbol,1],na.last=TRUE)), ]
write.table(tosave.m, file="TABLES/summary_mouse_allOrthologs.csv", row.names=F, quote=F, sep=",")

tosave.h <- formatResults(res.top.h)
tosave.h <- tosave.h[names(sort(res.top.h$Breakpoints[ortho.genes.use$hgnc_symbol,1],na.last=TRUE)), ]
write.table(tosave.h, file="TABLES/summary_human_allOrthologs.csv", row.names=F, quote=F, sep=",")



# Which genes start trend is either UP or DOWN and get INITIAL TIME that it starts.
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) which(x != 0)[1])]

sum(res.top.h$Trends[ortho.genes.use$hgnc_symbol,1] !=0) / length(ortho.genes.use$hgnc_symbol) # 90%
sum(res.top.m$Trends[ortho.genes.use$mgi_symbol,1] !=0) / length(ortho.genes.use$mgi_symbol) # 88%

X <- timeUpDown.h
Y <- timeUpDown.m

library(ggplot2)

pdf("PLOTS/histogram_FirstTime_UpDown_Figure2_Orthologs.pdf", height=6, width=8)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,100), 
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=50), 
	main="", xlab="Minute",
	cex.axis=2, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), breaks = seq(0, 600, length.out=50))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2, bty='n')
dev.off()


## Now get what percent of genes go up/down at time 0 versus delayed up/down
timeUpDown.m <- t.v.m[apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) which(x != 0)[1])]
timeUpDown.h <- t.v.h[apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) which(x != 0)[1])]

pcntStart0.h <- mean(timeUpDown.h==0)*100
pcntStart0.m <- mean(timeUpDown.m==0)*100

prop.test(c(sum(timeUpDown.h==0), sum(timeUpDown.m==0)), c(length(timeUpDown.h), length(timeUpDown.m)))
prop.test(c(sum(timeUpDown.h!=0), sum(timeUpDown.m!=0)), c(length(timeUpDown.h), length(timeUpDown.m)))

pdf("PLOTS/percent_FirstTime_UpDown_allTrendyGenes_Figure2_Orthologs.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(pcntStart0.m, pcntStart0.h, 100-pcntStart0.m, 100-pcntStart0.h),
space=c(.5,.1,1,.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,100), cex.axis=2
)
legend('topright', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()




## ORTHOLOGS Monotonic genes, only up or down the entire time:

onlyDown.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == -1)))
onlyDown.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == -1)))

onlyUp.m <- sum(apply(res.top.m$Trends[ortho.genes.use$mgi_symbol,], 1, function(x) all(x == 1)))
onlyUp.h <- sum(apply(res.top.h$Trends[ortho.genes.use$hgnc_symbol,], 1, function(x) all(x == 1)))

prop.test(c((onlyDown.m), (onlyDown.h)), c(length(ortho.genes.use$mgi_symbol), length(ortho.genes.use$mgi_symbol)))
prop.test(c((onlyUp.m), (onlyUp.h)), c(length(ortho.genes.use$mgi_symbol), length(ortho.genes.use$mgi_symbol)))

all.mono.m <- onlyDown.m + onlyUp.m
all.mono.h <- onlyDown.h + onlyUp.h

prop.test(c((all.mono.m), (all.mono.h)), c(length(ortho.genes.use$mgi_symbol), length(ortho.genes.use$mgi_symbol)))

pdf("PLOTS/numGenes_Monotonic_allTrendyGenes_Figure2_Orthologs.pdf", height=6, width=5)
par(mar=c(3,3,3,1))
barplot(c(all.mono.m, all.mono.h),
space=c(.1), col = c("cornflowerblue", "brown1"), names="",
ylim=c(0,50), cex.axis=2
)
legend('topleft', c("Mouse","Human"), fill=c("cornflowerblue", "brown1"), cex=2, bty='n')
dev.off()





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

X <- (na.omit(c(all.bp.human[,1])))
Y <- (na.omit(c(all.bp.mouse[,1])))

PP <- round(wilcox.test(X,Y)$p.value, 3)
if( PP < .001) {PP <- "< .001"}
PP
sum(X <= 250) / length(X)
sum(Y <= 250) / length(Y)
  
pdf("PLOTS/histogram_firstBreakpoint_Orthologs_Figure4.pdf", height=8, width=12)
par(mar=c(6,6,3,1), mgp=c(4,1,0))
hist(X, xlim=c(0,600), ylim=c(0,15), border="brown3",
	col=alpha("brown1", .6), breaks = seq(0, 600, length.out=20), 
	main="", xlab="Minute",
	cex.axis=2.5, cex.lab=3)
hist(Y, add=T,  col=alpha("cornflowerblue", .6), border="dodgerblue3", breaks = seq(0, 600, length.out=20))
legend('topright', c("Mouse","Human"), fill=c(alpha("cornflowerblue", .6), alpha("brown1", .6)), cex=2)
dev.off()

pdf("PLOTS/spectrum_firstBreakpoint_Orthologs_Figure4.pdf", height=3.5, width=16)
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

pdf("PLOTS/boxPlot_firstBreakpoint_Orthologs_Figure4.pdf", height=7, width=5)
par(mar=c(5,6,2,1), mgp = c(4, .5, 0))
pirateplot(formula = Minute ~ Species,
           data = longdata,
           xlab = "", inf.b.o = .3,point.o = .5,
           ylab = "Minutes", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=2.5, cex.axis=2,cex.names=2.5)
dev.off()


# mean(res.top.m$Breakpoints[ortho.genes.use$mgi_symbol,1], na.rm=T)
# mean(res.top.h$Breakpoints[ortho.genes.use$hgnc_symbol,1], na.rm=T)

mean(res.top.m$Breakpoints[ortho.genes.use$mgi_symbol,1] - res.top.h$Breakpoints[ortho.genes.use$hgnc_symbol,1], na.rm=T)



########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
