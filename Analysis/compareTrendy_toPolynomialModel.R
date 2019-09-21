setwd("~/RobotSeq/")

# Mouse original and scaled data:
load("RDATA/normalized_mouseRobot.RDATA")
data.norm.filter <- data.norm.mouse[which(rowMeans(data.norm.mouse) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
data.norm.scale.m <- data.norm.scale


load("RDATA/trendy_run_mouse_Scaled0to1.RData")
library(Trendy)
seg.mouse <- results(seg.all.scaled)
t.v.m <- t.v 


# # Get ALL genes fitted
res.top.m <- topTrendy(seg.mouse, -Inf)
length(res.top.m$AdjustedR2)

checkG <- names(res.top.m$AdjustedR2)

length(checkG)

fitall.r2 <- sapply(1:length(checkG), function(j){
  tt = lm(data.norm.scale.m[checkG[j],] ~ poly(t.v.m, 5))
  t2 = summary(tt)$adj.r.squared
  return(t2)  
})
names(fitall.r2) <- checkG

topSet <- names(which(res.top.m$AdjustedR2 < .2))
bottomSet <- names(which(res.top.m$AdjustedR2 > .2))
library(ggplot2)

pdf("PLOTS/forReviewer/hist_PolynomialR2_less.2-mouse.pdf",height=5, width=6)
par(mfrow=c(2,2), mar=c(5,5,2,1))
hist(res.top.m$AdjustedR2[topSet],  xlim=c(-.05,.6), xlab=expression(paste("Adjusted R"^"2")), main= "Mouse: Trendy Fit")
abline(v=.2, lwd=2, col="blue")
hist(fitall.r2[topSet], xlim=c(-.05,.6), xlab=expression(paste("Adjusted R"^"2")), main= "Mouse: Polynomial Fit")
abline(v=.2, lwd=2, col="blue")
polygon(c(.2, .2, .6, .6), c(0, 1500, 1500, 0),
     col = alpha("blue", .2), border = NA)
		 
hist(res.top.m$AdjustedR2[bottomSet],  xlim=c(-.05,1), xlab=expression(paste("Adjusted R"^"2")), main= "Mouse: Trendy Fit")
abline(v=.2, lwd=2, col="blue")
hist(fitall.r2[bottomSet], xlim=c(-.05,1), xlab=expression(paste("Adjusted R"^"2")), main= "Mouse: Polynomial Fit")
abline(v=.2, lwd=2, col="blue")
polygon(c(0, 0, .2, .2), c(0, 1500, 1500, 0),
     col = alpha("blue", .2), border = NA)
dev.off()

sum(fitall.r2[topSet] > .2)
sum(fitall.r2[bottomSet] < .2)



#########################################################



# Mouse original and scaled data:
load("RDATA/normalized_humanRobot.RDATA")
data.norm.filter <- data.norm.human[which(rowMeans(data.norm.human) >= 10), ]
data.norm.scale <- t(apply(data.norm.filter, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
data.norm.scale.h <- data.norm.scale


load("RDATA/trendy_run_human_Scaled0to1.RData")
library(Trendy)
seg.human <- results(seg.all.scaled)
t.v.h <- t.v 

# # Get ALL genes fitted
res.top.h <- topTrendy(seg.human, -Inf)
length(res.top.h$AdjustedR2)

checkG <- names(res.top.h$AdjustedR2)

length(checkG)

fitall.r2 <- sapply(1:length(checkG), function(j){
  tt = lm(data.norm.scale.h[checkG[j],] ~ poly(t.v.h, 5))
  t2 = summary(tt)$adj.r.squared
  return(t2)  
})
names(fitall.r2) <- checkG

topSet <- names(which(res.top.h$AdjustedR2 < .2))
bottomSet <- names(which(res.top.h$AdjustedR2 > .2))
library(ggplot2)

pdf("PLOTS/forReviewer/hist_PolynomialR2_less.2-human.pdf",height=5, width=6)
par(mfrow=c(2,2), mar=c(5,5,2,1))
hist(res.top.h$AdjustedR2[topSet],  xlim=c(-.05,.6), xlab=expression(paste("Adjusted R"^"2")), main= "Human: Trendy Fit")
abline(v=.2, lwd=2, col="blue")
hist(fitall.r2[topSet], xlim=c(-.05,.6), xlab=expression(paste("Adjusted R"^"2")), main= "Human: Polynomial Fit")
abline(v=.2, lwd=2, col="blue")
polygon(c(.2, .2, .6, .6), c(0, 1500, 1500, 0),
     col = alpha("blue", .2), border = NA)
		 
hist(res.top.h$AdjustedR2[bottomSet],  xlim=c(-.05,1), xlab=expression(paste("Adjusted R"^"2")), main= "Human: Trendy Fit")
abline(v=.2, lwd=2, col="blue")
hist(fitall.r2[bottomSet], xlim=c(-.05,1), xlab=expression(paste("Adjusted R"^"2")), main= "Human: Polynomial Fit")
abline(v=.2, lwd=2, col="blue")
polygon(c(0, 0, .2, .2), c(0, 1500, 1500, 0),
     col = alpha("blue", .2), border = NA)
dev.off()

sum(fitall.r2[topSet] > .2)
sum(fitall.r2[bottomSet] < .2)
