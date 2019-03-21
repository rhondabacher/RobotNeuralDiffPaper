setwd("~/RobotSeq/")


# Give objects distinct names:
load("RDATA/trendy_run_mouse_Scaled0to1.RData")

seg.mouse <- seg.all.scaled
data.norm.m <- data.norm.mouse
data.norm.scale.m <- data.norm.scale
t.v.m <- t.v 

library(Trendy)
peak_genes <- extractPattern(seg.mouse, Pattern=c("up", "down"), adjR2Cut =.5)
peak_genes.m <- peak_genes

load("RDATA/trendy_run_human_Scaled0to1.RData")

seg.human <- seg.all.scaled
data.norm.h <- data.norm.human
data.norm.scale.h <- data.norm.scale
t.v.h <- t.v 

peak_genes <- extractPattern(seg.human, Pattern=c("up", "down"), adjR2Cut =.5)
peak_genes.h <- peak_genes


## Get list of orthologs for comparing Mouse and Human:
library(biomaRt)


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes = c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_perc_id_r1")
attributes = c(attributes,"mmusculus_homolog_orthology_type", "mmusculus_homolog_subtype", "mmusculus_homolog_perc_id")

orth.mouse = getBM(attributes,filters="with_mmusculus_homolog",values =TRUE, mart = human, bmHeader=FALSE)

orth.mouse.unique <- orth.mouse
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human_convert <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", orth.mouse.unique$ensembl_gene_id, mart= human)
mouse_convert <- getBM(c("mgi_symbol","ensembl_gene_id"), "ensembl_gene_id", orth.mouse.unique$mmusculus_homolog_ensembl_gene, mart= mouse)

orth.1 <- merge(orth.mouse.unique, human_convert, by="ensembl_gene_id", all=TRUE)
orth.2 <- merge(orth.1, mouse_convert, by.x="mmusculus_homolog_ensembl_gene", by.y="ensembl_gene_id", all.x=TRUE)


orth.genes <- orth.2[,7:8]

# Clean this up a bit for future use:
orth.genes.clean1 <- orth.genes[which(orth.genes[,1]!="" & orth.genes[,2]!=""),]
orth.genes.clean <- orth.genes.clean1[!duplicated(orth.genes.clean1),]

# For using a cutoff of .5:
res.top.m <- topTrendy(seg.mouse, .5)
res.top.h <- topTrendy(seg.human, .5)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

top1 <- merge(orth.genes.clean, top.human, by="hgnc_symbol")
top2 <- merge(top1,top.mouse, by="mgi_symbol")
top2 <- top2[!duplicated(top2),]
dupg <- top2[which(duplicated(top2[,2])),2] 
subset(top2, hgnc_symbol %in% dupg)
dupg <- top2[which(duplicated(top2[,1])),1] 
subset(top2, mgi_symbol %in% dupg)
TORM <- c(45, 32)
top2 <- top2[-TORM,]
dim(top2)

# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse, !(Gene %in% top2$mgi_symbol))
intersect(toupper(mouse.genes.check1[,1]), top.human$Gene)
recovered.g <- data.frame(mgi_symbol = c("Ubc", "Iffo2"), hgnc_symbol = c("UBC","IFFO2"), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)

ortho.genes.use <- top2



save.image("RDATA/jointPlots_loadDataBoth.Rdata")