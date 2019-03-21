## All code relevant to Figure 1 plots and numbers:


setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")
library(Trendy)

# Mouse Plots

res.top <- topTrendy(seg.mouse, .5)

# How many breakpoints happen before 30 minutes?
sum(res.top$Breakpoint <= 30, na.rm=T)
# 83
sum(res.top$Breakpoint <= 100, na.rm=T) #264
sum(res.top$Breakpoint <= 250, na.rm=T) #738


# Heatmap of all genes with an adjusted R^2 >= .5
library(gplots)
heatcol <- colorpanel(100, "darkblue", "white", "darkred")
heatcol2 <- heatcol[c(1:20,seq(21,80,2),81:100)]
colnames(data.norm.scale) <- t.v

topg <- names(res.top$AdjustedR2)
MAX.H <- apply(data.norm.scale.m[topg,], 1, which.max)
data.norm.scale.sort <- data.norm.scale.m[names(sort(MAX.H)),]
colnames(data.norm.scale.sort) <- t.v

pdf("PLOTS/Mouse_ExprHeatMap_OrderByFirstBreakpoint.pdf", height=10)
par(mar=c(1,1,1,5), cex.axis=.5, cex.main=.5, cex.lab=.5, cex.sub=.5)
  heatmap.2(data.norm.scale.sort,
            trace="none", 
            Rowv=F,Colv=F,col=heatcol2, cex.main=.5, cex.axis=.5,
            cexRow=.1, cexCol = .8, key=TRUE, keysize=.5, lwid=c(1,4), lhei=c(1,4))
dev.off()


## Table for supplement:
tosave <- formatResults(res.top)
tosave <- tosave[rownames(data.norm.scale.sort), ] #order by time of first breakpoint
write.table(tosave, file="TABLES/Mouse_Fig1HeatmapGenes.csv", row.names=F, quote=F, sep=",")




##################################################
##################################################



# Human

# All genes with an adjusted R^2 >= .5
res.top <- topTrendy(seg.human, .5)

# How many breakpoints happen before 30 minutes?
sum(res.top$Breakpoint <= 30, na.rm=T)
# 0
sum(res.top$Breakpoint <= 100, na.rm=T) # 38
sum(res.top$Breakpoint <= 250, na.rm=T) # 100

# Heatmap of all genes with an adjusted R^2 >= .5
library(gplots)
heatcol <- colorpanel(100, "darkblue", "white", "darkred")
heatcol2 <- heatcol[c(1:20,seq(21,80,2),81:100)]
colnames(data.norm.scale) <- t.v

topg <- names(res.top$AdjustedR2)
MAX.H <- apply(data.norm.scale.h[topg,], 1, which.max)
data.norm.scale.sort <- data.norm.scale.h[names(sort(MAX.H)),]
colnames(data.norm.scale.sort) <- t.v

pdf("PLOTS/Human_ExprHeatMap_OrderByFirstBreakpoint.pdf", height=10)
par(mar=c(1,1,1,5), cex.axis=.5, cex.main=.5, cex.lab=.5, cex.sub=.5)
  heatmap.2(data.norm.scale.sort,
            trace="none", 
            Rowv=F,Colv=F,col=heatcol2, cex.main=.5, cex.axis=.5, 
            cexRow=.1, cexCol = .8, key=TRUE, keysize=.5, 
						density.info=c("none"),
						lwid=c(1,4), lhei=c(1,4))
dev.off()

## Table for supplement:
tosave <- formatResults(res.top)
tosave <- tosave[rownames(data.norm.scale.sort), ] #order by time of first breakpoint
write.table(tosave, file="TABLES/Human_Fig1HeatmapGenes.csv", row.names=F, quote=F, sep=",")




###################################################################################################
###################################################################################################
###################################################################################################

## Joint plots

# Dynamic genes in common:
library(Trendy)
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

sum(res.top.m$Breakpoint <= 30, na.rm=T)
sum(res.top.m$Breakpoint <= 100, na.rm=T)
sum(res.top.m$Breakpoint <= 250, na.rm=T)

sum(res.top.h$Breakpoint <= 30, na.rm=T)
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

pdf("PLOTS/boxPlot_GenesInFig1_AllBreakpoints.pdf", height=8, width=6)
par(mar=c(4,6,2,1), mgp = c(4, 1, 0))
pirateplot(formula = Breakpoint ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "Minute", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=3, cex.axis=2,cex.names=3)
dev.off()

X <- (na.omit(c(res.top.h$Breakpoints)))
Y <- (na.omit(c(res.top.m$Breakpoints)))


pdf("PLOTS/histogram_GenesInFig1_AllBreakpoints.pdf", height=12, width=15)
par(mar=c(5,5,4,2), mfrow=c(2,1))
hist(Y, xlim=c(0,600), ylim=c(0,50), 
	col="cornflowerblue", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="", xlab="",border="dodgerblue3", 
	cex.axis=2, cex.lab=3)
hist(X, xlim=c(0,600), ylim=c(0,50), 
	col="brown1", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="", xlab="",border="brown3", 
	cex.axis=2, cex.lab=3)
dev.off()



