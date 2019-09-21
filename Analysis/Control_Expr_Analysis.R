setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")
load("RDATA/trendy_run_control_Scaled0to1.Rdata")

## This code also generates Supplementary Figures XX

library(Trendy)

seg.human.control <- results(seg.all.scaled.control)

## Genes in common with control
top.res.h <- topTrendy(seg.human, adjR2Cut = .2)
top.res.h.control <- topTrendy(seg.human.control, adjR2Cut = .2)

str(top.res.h) # 4332
str(top.res.h.control) # 3327

genes.in.common <- intersect(names(top.res.h$AdjustedR2), names(top.res.h.control$AdjustedR2))
length(genes.in.common)

X <- (na.omit(top.res.h$Breakpoints[genes.in.common,1]))
Y <- (na.omit(top.res.h.control$Breakpoints[genes.in.common,1]))

pdf("PLOTS/histogram_GenesInCommon_RobotControlHuman-Breakpoints.pdf", height=4, width=7, useDingbats=F)
par(mar=c(2,5,1,2), mfrow=c(2,1))
hist(X, xlim=c(0,600), ylim=c(0,150),
	col="cornflowerblue", breaks = seq(0, 600, length.out=100), 
	ylab="Number of breakpoints", main="", xlab="",border="dodgerblue3", yaxt='n',
	cex.axis=1, cex.lab=1)
axis(2, seq(0,150, by=50))
abline(v=c(60, 80,280,300,440,490), col="gray45", lty=2, lwd=2)

hist(Y, xlim=c(0,600), ylim=c(0,150),
	col="brown1", breaks = seq(0, 600, length.out=100),  yaxt='n',
	ylab="Number of breakpoints", main="", xlab="",border="brown3", 
	cex.axis=1, cex.lab=1)
axis(2, seq(0,150, by=50))
abline(v=c(150,180,280,310,400,440), col="gray45", lty=2, lwd=2)
dev.off()


## Similarity based on start

diff.score <- (top.res.h.control$Segment.Trends[genes.in.common,1] + top.res.h$Segment.Trends[genes.in.common,1])
same.dir <- names(which(abs(diff.score) == 2)) # exact same direction
diff.dir <- setdiff(names(which(abs(diff.score) == 0)), names(which(top.res.h.control$Segment.Trends[genes.in.common,1] == 0)))
others <- setdiff(genes.in.common, c(diff.dir,same.dir))


write.table(same.dir, file="TABLES/ControlExpr_SameDirectionGenes.txt", quote=F, row.names=F, col.names=F)
write.table(diff.dir, file="TABLES/ControlExpr_DiffDirectionGenes.txt", quote=F, row.names=F, col.names=F)


sub1 <- names(which(top.res.h.control$Breakpoints[genes.in.common,1] > 150 & 
  top.res.h.control$Breakpoints[genes.in.common,1] < 180))
write.table(sub1, file="TABLES/firstSet_Control.txt", quote=F, row.names=F, col.names=F)

sub2 <- names(which(top.res.h.control$Breakpoints[genes.in.common,1] > 280 & 
  top.res.h.control$Breakpoints[genes.in.common,1] < 310))
write.table(sub2, file="TABLES/secondSet_Control.txt", quote=F, row.names=F, col.names=F)

sub3 <- names(which(top.res.h.control$Breakpoints[genes.in.common,1] > 400 & 
  top.res.h.control$Breakpoints[genes.in.common,1] < 440))
write.table(sub3, file="TABLES/thirdSet_Control.txt", quote=F, row.names=F, col.names=F)
  
sub4 <- names(which(top.res.h$Breakpoints[genes.in.common,1] > 60 & 
  top.res.h$Breakpoints[genes.in.common,1] < 80))
write.table(sub4, file="TABLES/firstSet_Human.txt", quote=F, row.names=F, col.names=F)

sub5 <- names(which(top.res.h$Breakpoints[genes.in.common,1] > 280 & 
  top.res.h$Breakpoints[genes.in.common,1] < 300))
write.table(sub5, file="TABLES/secondSet_Human.txt", quote=F, row.names=F, col.names=F)

sub6 <- names(which(top.res.h$Breakpoints[genes.in.common,1] > 440 & 
  top.res.h$Breakpoints[genes.in.common,1] < 490))
write.table(sub6, file="TABLES/thirdSet_Human.txt", quote=F, row.names=F, col.names=F)

c(length(intersect(sub4, c(sub1,sub2,sub3))) / length(sub4),
length(intersect(sub5, c(sub1,sub2,sub3))) / length(sub5),
length(intersect(sub6, c(sub1,sub2,sub3))) / length(sub6))



length(same.dir) / length(diff.score)
length(diff.dir) / length(diff.score)
length(others)
length(diff.dir)
length(same.dir)

library(ggsunburst)

df <- data.frame(
  parent = c("Trendy Genes", "Trendy Genes", "Shared", "Shared", "Shared"),
	node = c("Unique", "Shared", "Similar Direction", "Opposite Direction", "Others"),
  size = round(c((4332-1634)/4332, 1634/4332,  436/4332, 958/4332, 240/4332)*100,1)
  )

# write data.frame into csv file
write.table(df, file = 'df.csv', row.names = F, sep = ",")

# generate data structure
sb <- sunburst_data('df.csv', type = 'node_parent', sep = ",", node_attributes = c("size"))


pdf("PLOTS/pieChart_sharedGenes_SuppFig2.pdf", height=10, width=10, useDingbats=F)
par(mar=c(1,2,2,1))
p <- sunburst(sb, rects.fill.aes = "name", node_labels = T) +  scale_fill_brewer(palette = "Set3")  
p + geom_text(data = sb$leaf_labels,
    aes(x=x, y=0.1, label=paste(size,"%"), angle=0, hjust=1, vjust=-.5), size = 5)+ 
		geom_text(data = sb$node_labels,
    aes(x=x, y=0.1, label=label, vjust=c(1,-4)), size = 5)+ guides(fill=FALSE)
dev.off()
