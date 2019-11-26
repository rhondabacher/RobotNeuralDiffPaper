setwd("~/RobotNeuralDiffPaper/")


####### Mouse #########

load("RDATA/normalized_mouseRobot.RDATA")

cnames <- colnames(data.norm.mouse)
h.ind <- as.numeric(sapply(strsplit(cnames, "_"), function(x) x[2]))

t.v <- h.ind
names(t.v) <- cnames

############# Run Trendy on Raw Data ##########
library(Trendy)
BiocParallel::register(BiocParallel::SerialParam())

seg.all.orig <- trendy(data.norm.mouse, maxK = 5, tVectIn=t.v, pvalCut = .2, minNumInSeg=5,
   saveObject=TRUE, fileName="OUT/mouseRobot", numTry=5)

save(seg.all.orig, t.v, file="RDATA/trendy_run_mouse.RData")

############# Run Trendy on the Scaled Data ##########

# Remove lowly expressed genes as Trendy does by default:
data.norm.filter <- data.norm.mouse[which(rowMeans(data.norm.mouse) >= 10), ]

# Scale between 0 and 1:
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
BiocParallel::register(BiocParallel::SerialParam())

# Adjust the mean filter since we have scaled data:
seg.all.scaled <- trendy(data.norm.scale, meanCut=-Inf,  maxK = 5, tVectIn=t.v, pvalCut = .2, minNumInSeg=5,
   saveObject=TRUE, fileName="OUT/mouseRobot_scaled0to1", numTry=5)

save(seg.all.scaled, t.v, file="RDATA/trendy_run_mouse_Scaled0to1.RData")


############################################################################################################


####### Human #########

load("RDATA/normalized_humanRobot.RDATA")

cnames <- colnames(data.norm.human)
h.ind <- as.numeric(sapply(strsplit(cnames, "_"), function(x) x[2]))

t.v <- h.ind
names(t.v) <- cnames

############# Run Trendy ##########
library(Trendy)
BiocParallel::register(BiocParallel::SerialParam())

seg.all.orig <- trendy(data.norm.human, maxK = 5, tVectIn=t.v, pvalCut = .2, minNumInSeg=3,
   saveObject=TRUE, fileName="OUT/humanRobot", numTry=5)

save(seg.all.orig, t.v, file="RDATA/trendy_run_human.RData")

############# Run Trendy on the Scaled Data ##########

# Remove lowly expressed genes as Trendy does by default:
data.norm.filter <- data.norm.human[which(rowMeans(data.norm.human) >= 10), ]

# Scale between 0 and 1:
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
BiocParallel::register(BiocParallel::SerialParam())

# Adjust the mean filter since we have scaled data:
seg.all.scaled <- trendy(data.norm.scale, meanCut=-Inf, maxK = 5, tVectIn=t.v, pvalCut = .2, minNumInSeg=3,
   saveObject=TRUE, fileName="OUT/humanRobot_scaled0to1", numTry=5)

save(seg.all.scaled, t.v, file="RDATA/trendy_run_human_Scaled0to1.RData")



