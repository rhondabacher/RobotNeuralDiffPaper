setwd("~/RobotNeuralDiffPaper/")
library(Trendy)

# Give objects distinct names for mouse and human diff experiments:
load("RDATA/normalized_mouseRobot.RDATA")
data.norm.scale.m <- t(apply(data.norm.mouse, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

load("RDATA/trendy_run_mouse_Scaled0to1.RData")
seg.mouse <- results(seg.all.scaled)
t.v.m <- t.v 

load("RDATA/normalized_humanRobot.RDATA")
data.norm.scale.h <- t(apply(data.norm.human, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

load("RDATA/trendy_run_human_Scaled0to1.RData")
seg.human <- results(seg.all.scaled)
t.v.h <- t.v 

# Get genes that have ANY peak type pattern
peak_genes <- extractPattern(seg.mouse, Pattern=c("up", "down"), adjR2Cut = .2)
peak_genes <- rbind(peak_genes, extractPattern(seg.mouse, Pattern=c("up", "same", "down"), adjR2Cut = .2)[,1:2])
peak_genes.m <- peak_genes

peak_genes <- extractPattern(seg.human, Pattern=c("up", "down"), adjR2Cut = .2)
peak_genes <- rbind(peak_genes, extractPattern(seg.human, Pattern=c("up", "same", "down"), adjR2Cut = .2)[,1:2])
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

orth.genes <- orth.2[,6:8]

# Clean this up a bit for future use:
orth.genes.clean1 <- orth.genes[which(orth.genes[,1]!="" & orth.genes[,2]!=""),]
orth.genes.clean <- orth.genes.clean1[!duplicated(orth.genes.clean1),]

# For using a cutoff of .2:
res.top.m <- topTrendy(seg.mouse, .2)
res.top.h <- topTrendy(seg.human, .2)

top.mouse <- data.frame(Gene=names(res.top.m$AdjustedR2), mgi_symbol=names(res.top.m$AdjustedR2), row.names = names(res.top.m$AdjustedR2), stringsAsFactors=FALSE)
top.human <- data.frame(Gene=names(res.top.h$AdjustedR2), hgnc_symbol=names(res.top.h$AdjustedR2), row.names = names(res.top.h$AdjustedR2), stringsAsFactors=FALSE)

top1 <- merge(orth.genes.clean, top.human, by="hgnc_symbol")
top2 <- merge(top1,top.mouse, by="mgi_symbol")
top2 <- top2[!duplicated(top2),][,1:3]

dupg <- top2[which(duplicated(top2[,2])),2] 
subset(top2, hgnc_symbol %in% dupg)
toRM <- c()
for(i in 1:length(dupg)) {

  tempSet <- subset(top2, hgnc_symbol %in% dupg[i])
  bestMatch <- which(tempSet[,2] == toupper(tempSet[,1]))[1]
  if (!is.na(bestMatch)) {
    toRM <- c(toRM, rownames(subset(top2, hgnc_symbol %in% dupg[i]))[-bestMatch])
  }
  else {
    maxOrth <- which.max(subset(top2, hgnc_symbol %in% dupg[i])[,3])
    toRM <- c(toRM, rownames(subset(top2, hgnc_symbol %in% dupg[i]))[-maxOrth])
  }
}
dupg <- top2[which(duplicated(top2[,1])),1] 
subset(top2, mgi_symbol %in% dupg)
for(i in 1:length(dupg)) {

  tempSet <- subset(top2, mgi_symbol %in% dupg[i])
  bestMatch <- which(tempSet[,2] == toupper(tempSet[,1]))[1]
  if (!is.na(bestMatch)) {
    toRM <- c(toRM, rownames(subset(top2, mgi_symbol %in% dupg[i]))[-bestMatch])
  }
  else {
    maxOrth <- which.max(subset(top2, mgi_symbol %in% dupg[i])[,3])
    toRM <- c(toRM, rownames(subset(top2, mgi_symbol %in% dupg[i]))[-maxOrth])
  }
}
toRM <- unique(toRM)

top2 <- top2[-as.numeric(toRM),]
dim(top2)


# Any others might be missing?
mouse.genes.check1 <- subset(top.mouse, !(Gene %in% top2$mgi_symbol))
mouse.genes.check1 <- mouse.genes.check1[which(toupper(mouse.genes.check1[,1]) %in% top.human$Gene),]
recovered.g <- data.frame(mgi_symbol = mouse.genes.check1[,2], hgnc_symbol = toupper(mouse.genes.check1[,2]), stringsAsFactors=F)
top2 <- top2[,1:2]
top2 <- rbind(top2, recovered.g)

ortho.genes.use <- top2

dim(ortho.genes.use)



save.image("RDATA/jointPlots_loadDataBoth.Rdata")