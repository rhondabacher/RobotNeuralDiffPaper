setwd("~/RobotSeq/")

## Normalized and process the mouse and human differentiation time-course datasets.

# Data downloaded from GSE129014:

## Mouse Data:
inData <- read.table("DATA/GSE129014_mm10genes.no_mt.ec.tab", stringsAsFactors=F, 
			header=T, row.names=1, sep="\t", comment.char="", quote="")
data.mat <- data.matrix(inData[,c(145, 1:144)])

# Normalize
library(EBSeq)
Sizes <- MedianNorm(data.mat)
data.norm.mouse <- GetNormalizedMat(data.mat, Sizes)

save(data.norm.mouse, file="RDATA/normalized_mouseRobot.RDATA")

# Write normalized counts tables to CSV.
write.csv(data.norm.mouse, file="OUT/mouse_normalizedCounts.csv", quote=F)



## Human Data:
inData <- read.table("DATA/GSE129014_hg19genes.no_mt.ec.tab", stringsAsFactors=F, 
			header=T, row.names=1, sep="\t", comment.char="", quote="")
dim(inData)
data.mat <- data.matrix(inData[,1:121])
# Add tech reps together. All time points have two except 450min.
poolData <- c()
k = 1
while(k <= 90) {
  poolData <- cbind(poolData, rowSums(inData[,k:(k+1)]))
  k = k + 2
  print(k)
}
poolData <- cbind(poolData, inData[,91])
k = 92
while(k <= 120) {
  poolData <- cbind(poolData, rowSums(inData[,k:(k+1)]))
  k = k + 2
  print(k)
}
colnames(poolData) <- paste0("NBplusNog_", seq(0, 600, by=10)) 

head(poolData)

library(EBSeq)
Sizes <- MedianNorm(poolData)
data.norm.human <- GetNormalizedMat(poolData, Sizes)

save(data.norm.human, file="RDATA/normalized_humanRobot.RDATA")


write.csv(data.norm.human, file="OUT/human_normalizedCounts.csv", quote=F)
