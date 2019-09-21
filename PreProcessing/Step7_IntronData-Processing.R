setwd("~/RobotSeq/")

inData <- read.table("DATA/genes.intron_hg38.tab", comment.char="", header=1, row.names=1, stringsAsFactors = F)
inData[1:5,1:5]

# Map the column names:
getcolNames <- read.table("DATA/GSE129014_hg19genes.no_mt.ec.tab", comment.char="", header=1, stringsAsFactors = F, sep="\t")
getcolNames[1:5,1:5]
dim(inData)
dim(getcolNames)
colnames(inData) <- colnames(getcolNames)[2:183]


inData.expr <- inData[,1:121]


poolData <- c()
k = 1
while(k <= 90) {
  print(colnames(inData.expr[,k:(k+1)]))
  poolData <- cbind(poolData, rowSums(inData.expr[,k:(k+1)]))
  k = k + 2
  print(k)
}
poolData <- cbind(poolData, inData.expr[,91])
k = 92
while(k <= 120) {
  print(colnames(inData.expr[,k:(k+1)]))
  poolData <- cbind(poolData, rowSums(inData.expr[,k:(k+1)]))
  k = k + 2
  print(k)
}

colnames(inData.expr)
newNames <- paste0("human_", unique(sapply(strsplit(colnames(inData.expr), "_", fixed=T), function(x) x[2])))

colnames(poolData) <- newNames

data.mat <- data.matrix(poolData)
library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)
write.csv(data.norm, file="OUT/human_intronReads_normalized.csv", quote=F)

save(data.norm, file="RDATA/normalized_humanRobot_reMap-Introns.RDATA")



dim(data.norm)
cnames <- colnames(data.norm)
h.ind <- as.numeric(gsub("min", "", sapply(strsplit(cnames, "_"), function(x) x[2])))

t.v <- h.ind
names(t.v) <- cnames;

############# segmented regression ##########
library(Trendy)
BiocParallel::register(BiocParallel::SerialParam())

seg.human.intron <- trendy(data.norm, tVectIn = t.v, maxK = 5, minNumInSeg = 3, pvalCut = .2, saveObject=TRUE, 
																									fileName="OUT/human_remap_intron", numTry = 5)

BiocParallel::register(BiocParallel::SerialParam())

save(seg.human.intron, t.v, file="RDATA/trendy_run_human_remap_intron.RData")

############# segmented regression on scale data first ##########

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

BiocParallel::register(BiocParallel::SerialParam())
seg.human.intron.scaled <- trendy(data.norm.scale, meanCut=-Inf, maxK = 5, tVectIn=t.v, minNumInSeg = 5, pvalCut = .2,
	 saveObject=TRUE, fileName="OUT/human_remap_intronScaled", numTry = 5)

data.norm.scale.human.intron <- data.norm.scale
t.v.human <- t.v
data.norm.human.intron <- data.norm
save(seg.human.intron.scaled, data.norm.scale.human.intron, 
		data.norm.human.intron, t.v.human, file="RDATA/trendy_run_human_remap_intronScaled.RData")





#######################################################################################################################################
#######################################################################################################################################


### Also for Mouse:
inData <- read.table("DATA/genes.intron_mm10.tab", comment.char="", header=1, row.names=1, stringsAsFactors = F)
inData[1:5,1:5]
getcolNames <- read.table("DATA/GSE129014_mm10genes.no_mt.ec.tab", comment.char="", header=1, stringsAsFactors = F, sep="\t")
getcolNames[1:5,1:5]
dim(inData)
dim(getcolNames)
colnames(inData) <- colnames(getcolNames)[2:146]

inData[1:5,1:5]
inData <- inData[,c(145,1:144)]

data.mat <- data.matrix(inData)
library(EBSeq)
sizes <- MedianNorm(data.mat)
data.norm <- GetNormalizedMat(data.mat, sizes)

write.csv(data.norm, file="OUT/mouse_intronReads_normalized.csv", quote=F)

save(data.norm, file="RDATA/normalized_mouseRobot_reMap-Introns.RDATA")



cnames <- colnames(data.norm)
h.ind <- as.numeric(gsub("min", "", sapply(strsplit(cnames, "_"), function(x) x[2])))

t.v <- h.ind
names(t.v) <- cnames;

############# segmented regress
ion ##########
library(Trendy)
BiocParallel::register(BiocParallel::SerialParam())

seg.mouse.intron <- trendy(data.norm, tVectIn = t.v, maxK = 5, minNumInSeg = 5, pvalCut = .2, saveObject=TRUE, fileName="OUT/mouse_remap_intron", numTry = 5)

BiocParallel::register(BiocParallel::SerialParam())

save(seg.mouse.intron, t.v, file="RDATA/trendy_run_mouse_remap_intron.RData")

############# segmented regression on scale data first ##########

data.norm.filter <- data.norm[which(rowMeans(data.norm) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

BiocParallel::register(BiocParallel::SerialParam())
seg.mouse.intron.scaled <- trendy(data.norm.scale, meanCut=-Inf, maxK = 5, tVectIn=t.v, minNumInSeg = 5, pvalCut = .2,
	 saveObject=TRUE, fileName="OUT/mouse_remap_intronScaled", numTry = 5)

data.norm.scale.mouse.intron <- data.norm.scale
t.v.mouse <- t.v
data.norm.mouse.intron <- data.norm
save(seg.mouse.intron.scaled, data.norm.scale.mouse.intron, 
	data.norm.mouse.intron, t.v.mouse, file="RDATA/trendy_run_mouse_remap_intronScaled.RData")
	
	
	
	




