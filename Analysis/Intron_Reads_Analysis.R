setwd("~/RobotSeq/")


load("RDATA/jointPlots_loadDataBoth.Rdata")

library(Trendy)
seg.mouse.exon <- seg.mouse
seg.human.exon <- seg.human
data.norm.scale.h.exon <- data.norm.scale.h
data.norm.scale.m.exon <- data.norm.scale.m

res.top.h.exon <- topTrendy(seg.human.exon, .2)
res.top.m.exon <- topTrendy(seg.mouse.exon, .2)

ortho.genes.use.orig <- ortho.genes.use



## Now the Mouse INTRON Trendy results:
library(Trendy)

load("RDATA/trendy_run_mouse_remap_intronScaled.RDATA")
seg.mouse.intron.scaled <- results(seg.mouse.intron.scaled)

peak_genes <- extractPattern(seg.mouse.intron.scaled, Pattern=c("up", "down"), adjR2Cut =.2)
peak_genes <- rbind(peak_genes, extractPattern(seg.mouse.intron.scaled, Pattern=c("up", "same", "down"), adjR2Cut =.2)[,1:2])
peak_genes.m.intron <- peak_genes


## Now the Human INTRON Trendy results:
load("RDATA/trendy_run_human_remap_intronScaled.RDATA")
seg.human.intron.scaled <- results(seg.human.intron.scaled)

peak_genes <- extractPattern(seg.human.intron.scaled, Pattern=c("up", "down"), adjR2Cut =.2)
peak_genes <- rbind(peak_genes, extractPattern(seg.human.intron.scaled, Pattern=c("up", "same", "down"), adjR2Cut =.2)[,1:2])
peak_genes.h.intron <- peak_genes


## Now I will make a CLEAN list of all the orthologs in INTRON

# For using a cutoff of .2:
res.top.m.intron <- topTrendy(seg.mouse.intron.scaled, .2)
res.top.h.intron <- topTrendy(seg.human.intron.scaled, .2)

top.mouse.intron <- data.frame(Gene=names(res.top.m.intron$AdjustedR2), mgi_symbol=names(res.top.m.intron$AdjustedR2), row.names = names(res.top.m.intron$AdjustedR2), stringsAsFactors=FALSE)
top.human.intron <- data.frame(Gene=names(res.top.h.intron$AdjustedR2), hgnc_symbol=names(res.top.h.intron$AdjustedR2), row.names = names(res.top.h.intron$AdjustedR2), stringsAsFactors=FALSE)


top1 <- merge(orth.genes.clean, top.human.intron, by="hgnc_symbol")
top2 <- merge(top1,top.mouse.intron, by="mgi_symbol")
top2 <- top2[!duplicated(top2),][,1:3]


# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse.intron, !(Gene %in% top2$mgi_symbol))
mouse.genes.check1 <- mouse.genes.check1[which(toupper(mouse.genes.check1[,1]) %in% top.human.intron$Gene),]
recovered.g <- data.frame(mgi_symbol = mouse.genes.check1[,2], hgnc_symbol = toupper(mouse.genes.check1[,2]), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)

ortho.genes.use <- top2

dim(ortho.genes.use)




############# # Above was pre-processing #####


## Analysis below:

library(Trendy)

###################
# How many are peak orthologs?

# Ortholog Peaks:
peak_genes.m.intron <- peak_genes.m.intron[which(!duplicated(peak_genes.m.intron[,1])),]
peak_genes.h.intron <- peak_genes.h.intron[which(!duplicated(peak_genes.h.intron[,1])),]

set1 <- subset(ortho.genes.use, hgnc_symbol %in% peak_genes.h.intron$Gene)
ortho.peaks <- subset(set1, mgi_symbol %in% peak_genes.m.intron$Gene)

dim(ortho.peaks)

write.table(ortho.peaks[,1], file="OUT/intronOrtholog_peaks.txt", quote=F, col.names=F, row.names=F)

###################

ortho.genes.use <- ortho.genes.use[(ortho.genes.use[,2] %in% names(seg.human.exon)),]
ortho.genes.use <- ortho.genes.use[(ortho.genes.use[,1] %in% names(seg.mouse.exon)),]

CHK1 <- ortho.genes.use[,1]
CHK2 <- ortho.genes.use[,2]


### Check Diff:
CHK <- CHK2
diff.intron <- rowMeans(log(data.norm.scale.human.intron[CHK,52:61]+1)) - rowMeans(log(data.norm.scale.human.intron[CHK,1:10]+1))
diff.exon <- rowMeans(log(data.norm.scale.h.exon[CHK,52:61]+1)) - rowMeans(log(data.norm.scale.h.exon[CHK,1:10]+1))

pdf("PLOTS/intron_human_expressionChange.pdf", height=3, width=3, useDingbats=F)
par(mar=c(4,4,2,1), mgp=c(3,2,1))
plot(diff.exon, diff.intron, cex=1, xlab=bquote(Delta~ "Exon Expression"), ylab=bquote(Delta~ "Intron Expression"),
      cex.lab=1, cex.axis=1, pch=19, bty='n', main="Human", cex.main=1)
text(-.3,.35, label=paste0("R = ",round(cor(diff.exon, diff.intron), 2)), cex=1)
dev.off()

CHK <- CHK1
diff.intron <- rowMeans(log(data.norm.scale.mouse.intron[CHK,136:145]+1)) - rowMeans(log(data.norm.scale.mouse.intron[CHK,1:10]+1))
diff.exon <- rowMeans(log(data.norm.scale.m.exon[CHK,136:145]+1)) - rowMeans(log(data.norm.scale.m.exon[CHK,1:10]+1))

pdf("PLOTS/intron_mouse_expressionChange.pdf", height=3, width=3, useDingbats=F)
par(mar=c(4,4,2,1), mgp=c(3,2,1))
plot(diff.exon, diff.intron, cex=1, xlab=bquote(Delta~ "Exon Expression"), ylab=bquote(Delta~ "Intron Expression"),
      cex.lab=1, cex.axis=1, pch=19, bty='n', main="Mouse", cex.main=1)
text(-.3,.2, label=paste0("R = ",round(cor(diff.exon, diff.intron), 2)), cex=1)
dev.off()




###################


upgenes.m <- c(names(which(res.top.m.exon$Segment.Trends[,1] == 0 & res.top.m.exon$Segment.Trends[,2] == 1)), 
                names(which(res.top.m.exon$Segment.Trends[,1] == 1)))

upgenes.h <- c(names(which(res.top.h.exon$Segment.Trends[,1] == 0 & res.top.h.exon$Segment.Trends[,2] == 1)), 
                names(which(res.top.h.exon$Segment.Trends[,1] == 1)))



all.slopes.mouse <- intersect(upgenes.m, names(res.top.m.intron$AdjustedR2))
all.slopes.mouse <- res.top.m.intron$Segment.Slopes[all.slopes.mouse,1]
all.slopes.human <- intersect(upgenes.h, names(res.top.h.intron$AdjustedR2))
all.slopes.human <- res.top.h.intron$Segment.Slopes[all.slopes.human,1]


library(pairwiseCI)
Y <- data.frame( Slope = (all.slopes.human[all.slopes.human>0]), Species = "Human")
X <- data.frame( Slope = (all.slopes.mouse[all.slopes.mouse>0]), Species = "Mouse")

longdata <- rbind(Y, X)
propTestVals <- pairwiseCI(Slope ~ Species, data = longdata, conf.level = 0.99, method='Median.ratio')
propTestVals <- round(do.call(c,propTestVals$byout[[1]][2:3]), 3)

X <- data.frame( Slope = (all.slopes.human[all.slopes.human>0]), Species = "Human")
Y <- data.frame( Slope = (all.slopes.mouse[all.slopes.mouse>0]), Species = "Mouse")
longdata <- rbind(Y, X)

pdf("PLOTS/boxPlot_anyUpExon_slopeUpIntron.pdf", height=3, width=3)
par(mar=c(1,2.7,1,.1), mgp=c(2,.5,0))
pirateplot(formula = Slope ~ Species,
	           data = longdata, avg.line.fun =median, avg.line.lwd=.8,avg.line.o=1,
	           xlab = "", inf.b.o = .5, point.o = .5, bar.f.o = 0, bean.f.o = 1,
	           ylab = "Slope", pal=c("cornflowerblue", "brown1"),inf.method = "iqr",
	           main = "", point.cex=.3,cex.lab=.7, cex.axis=.6,cex.names=.7)
mtext(c("Mouse", "Human"), side=1, at = c(1,2), cex=.6)
mtext(bquote("99% CI " ~ Delta~ "S: ("~ .(propTestVals[1]) ~ ", "~ .(propTestVals[2]) ~ ")" ), side=3, at = c(1.5), cex=.6)
dev.off()
