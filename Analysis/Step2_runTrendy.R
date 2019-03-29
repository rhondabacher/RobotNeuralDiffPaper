setwd("~/RobotSeq/")


####### Mouse #########

load("RDATA/normalized_mouseRobot.RDATA")

cnames <- colnames(data.norm.mouse)
h.ind <- as.numeric(sapply(strsplit(cnames, "_"), function(x) x[2]))

t.v <- h.ind
names(t.v) <- cnames

############# Run Trendy ##########
library(Trendy)

seg.all.orig <- trendy(data.norm.mouse, Max.K = 5, T.Vect=t.v, Pval.Cut = .2, Save.Object=TRUE, File.Name="OUT/mouseRobot")

save.image("RDATA/trendy_run_mouse.RData")

############# Run Trendy on the Scaled Data ##########

# Remove lowly expressed genes as Trendy does by default:
data.norm.filter <- data.norm.mouse[which(rowMeans(data.norm.mouse) >= 10), ]

# Scale between 0 and 1:
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

# Adjust the mean filter since we have scaled data:
seg.all.scaled <- trendy(data.norm.scale, Mean.Cut=-Inf, Max.K = 5, T.Vect=t.v, Pval.Cut = .2, Save.Object=TRUE, File.Name="OUT/Mouse/mouseRobot_scaled0to1")

save.image("RDATA/trendy_run_mouse_Scaled0to1.RData")


############################################################################################################


####### Human #########

load("RDATA/normalized_humanRobot.RDATA")

cnames <- colnames(data.norm.human)
h.ind <- as.numeric(sapply(strsplit(cnames, "_"), function(x) x[2]))

t.v <- h.ind
names(t.v) <- cnames

############# Run Trendy ##########
library(Trendy)

seg.all.orig <- trendy(data.norm.human, Max.K = 5, T.Vect=t.v, Pval.Cut = .2, Save.Object=TRUE, File.Name="OUT/humanRobot")

save.image("RDATA/trendy_run_human.RData")

############# Run Trendy on the Scaled Data ##########

# Remove lowly expressed genes as Trendy does by default:
data.norm.filter <- data.norm.human[which(rowMeans(data.norm.human) >= 10), ]

# Scale between 0 and 1:
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

# Adjust the mean filter since we have scaled data:
seg.all.scaled <- trendy(data.norm.scale, Mean.Cut=-Inf, Max.K = 5, T.Vect=t.v, Pval.Cut = .2, Save.Object=TRUE, File.Name="OUT/humanRobot_scaled0to1")

save.image("RDATA/trendy_run_human_Scaled0to1.RData")



