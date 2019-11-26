setwd("~/RobotNeuralDiffPaper/")

load("RDATA/normalized_mouseRobot.RDATA")
dim(data.norm.mouse)
head(data.norm.mouse)

cnames <- colnames(data.norm.mouse)
h.ind <- as.numeric(sapply(strsplit(cnames, "_"), function(x) x[2]))

t.v <- h.ind
names(t.v) <- cnames;


## Subsample to every 12 hours.
data.norm.mouse.short <- data.norm.mouse[,seq(1,147, by = 3)]
t.v <- t.v[seq(1,147, by = 3)]


############# segmented regression on scale data  ##########
library(Trendy)
data.norm.filter <- data.norm.mouse.short[which(rowMeans(data.norm.mouse.short) >= 10), ]
data.norm.mouse.short.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

BiocParallel::register(BiocParallel::MulticoreParam())
seg.mouse.scaled.short <- trendy(data.norm.mouse.short.scale, meanCut=-Inf,
		 												maxK = 5, tVectIn=t.v, minNumInSeg = 5, pvalCut = .2, 
										saveObject=TRUE, fileName="OUT/mouseRobot_scaleData_shortSubset", numTry = 5)

t.v.short <- t.v
save(data.norm.mouse.short.scale, seg.mouse.scaled.short, t.v.short, file="RDATA/trendy_run_mouse_scaleData_shortSubset.RData")

#############
