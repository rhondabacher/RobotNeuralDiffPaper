setwd("~/RobotSeq/")


## Mouse Data:
inData <- read.csv("DATA/MOUSE_genes.no_mt.ec.csv", stringsAsFactors=F, header=T)

rownames(inData) <- inData[,1]
data.mat <- data.matrix(inData[,-1])

# Normalize
library(EBSeq)
Sizes <- MedianNorm(data.mat)
data.norm.mouse <- GetNormalizedMat(data.mat, Sizes)


save(data.norm.mouse, file="RDATA/normalized_mouseRobot.RDATA")


## Human Data:
inData <- read.csv("DATA/HUMAN_genes.no_mt.ec.csv", stringsAsFactors=F, header=T)

rownames(inData) <- inData[,1]
data.mat <- data.matrix(inData[,-1])


library(EBSeq)
Sizes <- MedianNorm(data.mat)
data.norm.human <- GetNormalizedMat(data.mat, sizes)

save(data.norm.human, file="RDATA/normalized_humanRobot.RDATA")


