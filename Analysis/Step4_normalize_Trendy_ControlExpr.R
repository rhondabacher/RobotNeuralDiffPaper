setwd("~/RobotSeq/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

library(readxl)
inFile <- "DATA/Sub_robot_hg19__6145aa5fa175a408.xlsx"
inData <- read_excel(inFile, sheet=3)

data.mat <- data.matrix(inData)[,c(123:183)]
Genes <- as.matrix(as.data.frame(inData[,1], stringsAsFactors=F))
rownames(data.mat) <- Genes
  
library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)

save(data.norm, file="RDATA/normalized_humanControl_Robot.RData")
head(data.norm)

cnames <- colnames(data.mat)
h.ind <- as.numeric(gsub("s","", gsub("min", "", sapply(strsplit(cnames, "_"), function(x) x[2]))))

t.v <- h.ind
names(t.v) <- cnames;
print(t.v)

############# Trendy ##########
library(Trendy)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all <- trendy(Data = data.norm.scale, tVectIn = t.v, saveObject = TRUE, fileName = "RDATA/humanControl_Scaled0to1", pvalCut = .2, maxK = 5, meanCut = -Inf)

ready.trendy <- results(seg.all)

save.image("RDATA/trendy_run_control_Scaled0to1.Rdata")


