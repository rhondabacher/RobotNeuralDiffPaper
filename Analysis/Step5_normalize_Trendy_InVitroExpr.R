setwd("~/RobotSeq/")

# Data from previously published Barry et al. 2017 in GEO: GSE90053

## Mouse Data:

library(readxl)
inFile <- "GSE90053_mEpiVH1.RSEM.xlsx"
inData <- read_excel(inFile, sheet=4)

data.mat <- data.matrix(inData)[,2:17]
Genes <- c(as.matrix(as.data.frame(inData[,1])))
rownames(data.mat) <- Genes

library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)

save(data.norm, file="RDATA/normalized_mouseInVitro.RData")

cnames <- colnames(data.mat)
cnames <- sapply(strsplit(cnames,"_"), function(x) rev(x)[1])
h.ind <- as.numeric(gsub("d", "", cnames))
  
t.v <- h.ind
names(t.v) <- cnames;
print(t.v)


############# Trendy ##########
library(Trendy)

seg.all <- trendy(Data = data.norm, tVectIn = t.v, saveObject = TRUE, fileName = "OUT/MouseInVitro", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = 10)
ready.trendy <- results(seg.all)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all.scaled <- trendy(Data = data.norm.scale, tVectIn = t.v, saveObject = TRUE, fileName = "OUT/MouseInVitro_scaled0to1", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = -Inf)

ready.trendy.scaled <- results(seg.all.scaled)

save.image("RDATA/trendy_run_mouse_InVitro.Rdata")


##############################################################################
##############################################################################
##############################################################################

## Human Data:

library(readxl)
inFile <- "GSE90053_mEpiVH1.RSEM.xlsx"
inData <- read_excel(inFile, sheet=2)

data.mat <- data.matrix(inData)[,2:27]
Genes <- c(as.matrix(as.data.frame(inData[,1])))
rownames(data.mat) <- Genes

library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)

save(data.norm, file="RDATA/normalized_humanInVitro.RData")

cnames <- colnames(data.mat)
h.ind <- as.numeric(substr(cnames,10,11))

t.v <- h.ind
names(t.v) <- cnames;
print(t.v)


############# Trendy ##########
library(Trendy)

seg.all <- trendy(Data = data.norm, tVectIn = t.v, saveObject = TRUE, fileName = "OUT/HumanInVitro", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = 10)
ready.trendy <- results(seg.all)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all.scaled <- trendy(Data = data.norm.scale, tVectIn = t.v, saveObject = TRUE, fileName = "OUT/HumanInVitro_scaled0to1", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = -Inf)

ready.trendy.scaled <- results(seg.all.scaled)

save.image("RDATA/trendy_run_human_InVitro.Rdata")


