setwd("~/RobotSeq/")


## Mouse Data:

library(readxl)
inFile <- "DATA/MouseInVitro.xlsx"
inData <- read_excel(inFile, sheet=1)

inData[which(duplicated(inData[,1])),1]
inData[which(inData[,1] == "42430"),1] <- c("42430.1", "42430.2")
inData[which(inData[,1] == "42431"),1] <- c("42431.1", "42431.2")

data.mat <- data.matrix(inData)[,2:17]
Genes <- c(as.matrix(as.data.frame(inData[,1])))

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

seg.all <- trendy(Data = data.norm, tVectIn = t.v, saveObject = TRUE, fileName = "RDATA/MouseInVitro", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = 10)
ready.trendy <- results(seg.all)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all.scaled <- trendy(Data = data.norm.scale, tVectIn = t.v, saveObject = TRUE, fileName = "RDATA/MouseInVitro_Scaled0to1", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = -Inf)

ready.trendy.scaled <- results(seg.all.scaled)

save.image("RDATA/trendy_run_mouse_InVitro.Rdata")


##############################################################################
##############################################################################
##############################################################################

## Human Data:
library(readxl)
inFile <- "DATA/HumanInVitro.xlsx"
inData <- read_excel(inFile, sheet=1)

inData[which(duplicated(inData[,1])),1]
inData[which(inData[,1] == "42431"),1] <- c("42431.1", "42431.2")
inData[which(inData[,1] == "42430"),1] <- c("42430.1", "42430.2")

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

seg.all <- trendy(Data = data.norm, tVectIn = t.v, saveObject = TRUE, fileName = "RDATA/HumanInVitro", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = 10)
ready.trendy <- results(seg.all)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all.scaled <- trendy(Data = data.norm.scale, tVectIn = t.v, saveObject = TRUE, fileName = "RDATA/HumanInVitro_Scaled", pvalCut = .2, maxK = 5, minNumInSeg = 3, meanCut = -Inf)

ready.trendy.scaled <- results(seg.all.scaled)

save.image("RDATA/trendy_run_human_InVitro.Rdata")


