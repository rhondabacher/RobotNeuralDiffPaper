setwd("~/RobotNeuralDiffPaper/")

load("RDATA/jointPlots_loadDataBoth.Rdata")



peaks.orth.time <- merge(ortho.genes.use, peak_genes.h, by.x="hgnc_symbol", by.y="Gene")
head(peaks.orth.time); dim(peaks.orth.time)
head(peak_genes.m)

peaks.orth.time <- merge(peaks.orth.time, peak_genes.m, by.x="mgi_symbol", by.y="Gene")
head(peaks.orth.time)
dim(peaks.orth.time)


# Plot the first peak time for each species:
peaks.orth.time <- peaks.orth.time[which(!duplicated(peaks.orth.time[,1])),]

pdf("PLOTS/peakTime_MouseVsHuman.pdf", height=3.5, width=4, useDingbats=FALSE)
par(mar=c(5,5,2,1))
plot(peaks.orth.time[,4], peaks.orth.time[,3], xlim=c(0,600), ylim=c(0,600),
xlab="Mouse Peak Time", ylab="Human Peak Time", pch=16, cex=.7, main="Initial Peak Times")
dev.off()



                                                            