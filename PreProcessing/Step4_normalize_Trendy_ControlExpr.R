setwd("~/RobotNeuralDiffPaper/")

load("RDATA/jointPlots_loadDataBoth.Rdata")

# Data downloaded from GSE129014:
library(readxl)
inData <- read.table("DATA/GSE129014_hg19genes.no_mt.ec.tab", stringsAsFactors=F, 
			header=T, row.names=1, sep="\t", comment.char="", quote="")
dim(inData)
data.mat <- data.matrix(inData[,122:182])
  
library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)

# Write normalized counts tables to CSV.
write.csv(data.norm, file="OUT/humanControl_normalizedCounts.csv", quote=F)


save(data.norm, file="RDATA/normalized_humanControl_Robot.RData")
head(data.norm)
data.norm.control <- data.norm

cnames <- colnames(data.mat)
h.ind <- as.numeric(gsub("s","", gsub("min", "", sapply(strsplit(cnames, "_"), function(x) x[2]))))

t.v <- h.ind
names(t.v) <- cnames;
print(t.v)

############# Trendy ##########
library(Trendy)

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

seg.all.scaled.control <- trendy(Data = data.norm.scale, tVectIn = t.v, minNumInSeg = 3, saveObject = TRUE, fileName = "OUT/humanRobotControl_Scaled0to1", 
				pvalCut = .2, maxK = 5, meanCut = -Inf, numTry=5)


data.norm.scale.control <- data.norm.scale

save(data.norm.control, data.norm.scale.control, seg.all.scaled.control, t.v, file="RDATA/trendy_run_control_Scaled0to1.Rdata")


