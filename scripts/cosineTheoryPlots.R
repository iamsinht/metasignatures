library(metasignatures)

outdir <- "~/Work/bhk/analysis/cosine/manuscriptFiles/"
figdir <- "~/Work/bhk/analysis/cosine/manuscriptFigures/"

#### Case 1 simulation 
if (!file.exists(file.path(outdir, "simulCase1CosVsDim.rds"))){
  cosDim <- cosineVsDim(dims=seq(1000))
  cosDim$theoryVar <- 1/cosDim$dim
  saveRDS(cosDim, file.path(outdir, "simulCase1CosVsDim.rds"))
} else {
  cosDim <- readRDS(file.path(outdir, "simulCase1CosVsDim.rds"))
}

pdf(file.path(figdir, "simulCase1CosVsDim.pdf"), width=8, height=6)
ggplot(cosDim, aes(x=sd^2, y=theoryVar, color=dim)) + geom_point() + theme_minimal() + 
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") + xlab("Empirical variance") + 
  ylab("Theoretical variance") + ggtitle("Case 1: Variance of cosine similarity for multivariate standard normals") + 
  geom_abline(slope = 1, intercept=0, lty=2, col="grey")
dev.off()


#### Approximation of normal length simulation

# Load the latest version from H4H; placeholder
if (!file.exists(file.path(outdir, "simulCase2NormLength_iter=3_distiter=10.rds"))){
  nLen <- getNormalLength(ndim=c(seq(100, 1000, 100), seq(1500, 3000, 500)), iter=3, distiter = 10)
  saveRDS(nLen, file.path(outdir, "simulCase2NormLength_iter=3_distiter=10.rds"))
} else {
  nLen <- readRDS(file.path(outdir, "simulCase2NormLength_iter=3_distiter=10.rds"))
}


pdf(file.path(figdir, "simulCase2MeanLengthConvergence.pdf"), width=8, height=6)
ggplot(nLen, aes(x=ndim, y=meanLen/theoryLength)) + geom_point(col="royalblue4") + theme_minimal() + xlab("Dimension") + 
  ylab("Mean Length / Theory") + ggtitle("Approximation of Length of multivariate normal vectors") + 
  ylim(c(0.9, 1.05)) + geom_hline(yintercept=1, lty=2, col="grey")
dev.off()

pdf(file.path(figdir, "simulCase2sdLenConvergence.pdf"), width=8, height=6)
ggplot(nLen, aes(x=ndim, y=sdLen/meanLen)) + geom_point(col="orangered3") + theme_minimal() + xlab("Dimension") +
  ylab("SD Length / Mean Length") + ggtitle("Convergence of SD Length/Mean Length for multivariate normal vectors")
dev.off()


#### Case 2 simulation



#### Case 3 simulation
if (!file.exists(file.path(outdir, "simulCase3CosVsDim.rds"))){
  gVar2 <- sampleGenVarCos(iter=1, distiter=100)
  saveRDS(gVar2, file.path(outdir, "simulCase3CosVsDim.rds"))
} else {
  gVar2 <- readRDS(file.path(outdir, "simulCase3CosVsDim.rds"))
}

pdf(file.path(figdir, "simulCase3CosVsDim.pdf"), width=8, height=6)
ggplot(gVar2, aes(x=-log10(obsVar), y=-log10(theoryVar), color=as.factor(ndim))) + geom_point() + 
  xlim(c(0, 12)) + ylim(c(0,12)) + geom_abline(intercept=0, slope=1, col="red", lty=2) + 
  theme_minimal() + xlab("-Log10 Observed Variance") + ylab("-Log10 Theoretical Variance") + 
  ggtitle(sprintf("Case 3 Variance simulation, Pearson = %0.4f", cor(gVar2$obsVar, gVar2$theoryVar)))
dev.off()