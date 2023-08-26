library(ggplot2)

plotMoments <- function(moment5dir, figdir="~/Work/bhk/analysis/metasig/figures"){ #figdir="."){
  
  # Level 5 moments
  LMmomentFiles <- list.files(moment5dir, pattern="LM.rds")
  
  allix <- grep("L1K_all_LM.rds", LMmomentFiles)
  
  allMoments <- readRDS(file.path(moment5dir, LMmomentFiles[allix]))
  allM10 <- readRDS(file.path(moment5dir, "metasigs/L1K_all_k=10_LM.rds"))
  allM100 <- readRDS(file.path(moment5dir, "metasigs/L1K_all_k=100_LM.rds"))
  
  
  pdf(file.path(figdir, "moments_skewVsMean_allL5_LM.pdf"), width=8, height=6)
  ggplot(allMoments, aes(x=skew, y=mean)) + geom_point() + theme_minimal() + xlab("Skewness") + ylab("Mean Z-score") + 
    ggtitle(sprintf("L1000 Landmark Mean vs Skewness, All data, Pearson=%0.4f", 
                    cor(allMoments$skew, allMoments$mean))) + xlim(c(-3.5, 3.5)) + ylim(c(-0.8, 0.8))
  dev.off()
  
  pdf(file.path(figdir, "moments_MeanVsSD_allL5_LM.pdf"), width=8, height=6)
  ggplot(allMoments, aes(x=mean, y=sqrt(var))) + geom_point() + theme_minimal() + xlab("Mean Z-score") + 
    ylab("SD Z-score") + ggtitle(sprintf("L1000 Landmark Mean vs Standard Dev, All data")) + 
    xlim(c(-0.8, 0.8)) + ylim(c(0, 2.5))
  dev.off()
  
  combdf <- rbind(cbind(allMoments, metasize=1), cbind(allM10, metasize=10), cbind(allM100, metasize=100))
  
  pdf(file.path(figdir, "moments_MeanVsSD_allL5_LM_Meta.pdf"), width=8, height=6)
  ggplot(combdf, aes(x=mean, y=sqrt(var), color=as.factor(metasize))) + geom_point() + theme_minimal() + xlab("Mean Z-score") + 
    ylab("SD Z-score") + ggtitle(sprintf("L1000 Landmark Gene Mean vs Standard Dev")) + 
    xlim(c(-0.8, 0.8)) + ylim(c(0, 2.5)) + geom_vline(xintercept=0, lty=2) + guides(color = guide_legend(title="Metasize"))
  dev.off()
  
  # Skew 95 is the skewness of the 5th - 95th percentiles of data, i.e. the inner 90 percent. 
  pdf(file.path(figdir, "moments_SkewVsSkew95_allL5_LM.pdf"), width=8, height=6)
  ggplot(allMoments, aes(x=skew, y=skew95)) + geom_point() + theme_minimal() + xlab("Skewness") + 
    ylab("Skewness of inner 90 percent of data") + geom_abline(intercept=0, slope=1, col="blue", lty=2) + 
    ggtitle(sprintf("L1000 Landmark gene skewness, Pearson = %0.3f", cor(allMoments$skew, allMoments$skew95))) + 
    xlim(c(-3, 3)) + ylim(c(-3,3))
  dev.off()
  
  
  totalMoments <- cbind(allMoments, dataset="all data")
  
  for (ii in setdiff(seq_along(LMmomentFiles), allix)){
    mycell = strsplit(LMmomentFiles[ii], "_")[[1]][2]
    
    cellMoments <- readRDS(file.path(moment5dir, LMmomentFiles[ii]))
    totalMoments <- rbind(totalMoments, cbind(cellMoments, dataset=mycell))
  } 
  
  pdf(file.path(figdir, "moments_skewVsMean_allCellLinesL5_LM.pdf"), width=8, height=6)
  ggplot(totalMoments, aes(x=skew, y=mean, color=dataset)) + geom_point(alpha=0.3) + theme_minimal() + xlab("Skewness") + ylab("Mean Z-score") + 
    ggtitle(sprintf("L1000 Landmark Mean vs Skewness, All cell lines, Pearson=%0.4f", 
                    cor(totalMoments$skew, totalMoments$mean))) + xlim(c(-4, 4)) + ylim(c(-1.2, 1.2))
  dev.off()
  
  pdf(file.path(figdir, "moments_MeanVsSD_allCellLinesL5_LM.pdf"), width=8, height=6)
  ggplot(totalMoments, aes(x=mean, y=sqrt(var), color=dataset)) + geom_point(alpha=0.3) + theme_minimal() + xlab("Mean Z-score") + 
    ylab("SD Z-score") + ggtitle(sprintf("L1000 Landmark Mean vs Standard Dev, All cell lines")) + 
    xlim(c(-1.65, 1.65)) + ylim(c(0, 3.5))
  dev.off()
  
  
  
}