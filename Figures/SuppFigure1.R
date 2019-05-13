setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")


load("RDATA/trendy_run_control_Scaled0to1.Rdata")
seg.human.control <- ready.trendy
data.norm.scale.h.control <- data.norm.scale
t.v.h.control <- t.v 


library(Trendy)

## Genes in common with control

top.res.h <- topTrendy(seg.human, adjR2Cut = .5)
top.res.h.control <- topTrendy(seg.human.control, adjR2Cut = .5)

genes.in.common <- intersect(names(top.res.h$AdjustedR2), names(top.res.h.control$AdjustedR2))
length(genes.in.common)

## Similarity based on start

diff.score <- (top.res.h.control$Segment.Trends[genes.in.common,1] + top.res.h$Segment.Trends[genes.in.common,1])

same.dir <- names(which(abs(diff.score) == 2))
diff.dir <- setdiff(names(which(abs(diff.score) == 0)), names(which(top.res.h.control$Segment.Trends[genes.in.common,1] == 0)))

others <- setdiff(genes.in.common, c(diff.dir,same.dir))



library(ggsunburst)

df <- data.frame(
  parent = c("Trendy Genes", "Trendy Genes", "Shared", "Shared", "Shared"),
	node = c("Unique", "Shared", "Similar Direction", "Opposite Direction", "Others"),
  size = round(c((655/739), (84/739),  29/739, 31/739, 24/739)*100,1)
  )

# write data.frame into csv file
write.table(df, file = 'df.csv', row.names = F, sep = ",")

# generate data structure
sb <- sunburst_data('df.csv', type = 'node_parent', sep = ",", node_attributes = c("size"))


pdf("PLOTS/pieChart_sharedGenes_SuppFig2.pdf", height=10, width=10)
par(mar=c(1,2,2,1))
p <- sunburst(sb, rects.fill.aes = "name", node_labels = T) +  scale_fill_brewer(palette = "Set3")  
p + geom_text(data = sb$leaf_labels,
    aes(x=x, y=0.1, label=paste(size,"%"), angle=0, hjust=1, vjust=-.5), size = 5)+ 
		geom_text(data = sb$node_labels,
    aes(x=x, y=0.1, label=label, vjust=c(1,-4)), size = 5)+ guides(fill=FALSE)
dev.off()