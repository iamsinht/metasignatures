library(ggplot2)
library(ggridges)
library(grid)

plotMetasigs <- function(figdir = "~/Work/bhk/analysis/metasig/figures", datadir = "~/Work/bhk/analysis/metasig/h4h"){
  
  #### Figure 1C: Self-correlation curves ####
  xbase <- readRDS(file.path(datadir, "base/L1KCP", "AllDSMetaSimpearson500x100.rds"))
  xdf <- dplyr::bind_rows(lapply(seq_len(length(xbase)), FUN=function(y) cbind(xbase[[y]]$cordf, compound=names(xbase)[y])))
  
  ndmsoBase <- readRDS(file.path(datadir, "base/dmso", "allDSDMSOMetaSimpearson500x100.rds"))
  nbgcpBase <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimPearson500x100.rds"))
  
  pdf(file.path(figdir, "L1KCPAll_n=100.pdf"), width=8, height=6)
  ggplot() + theme_minimal() +
    geom_line(data=xdf, aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoBase$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgcpBase$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound metacorrelation, unnormalized")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom", text=element_text(size=16))
  dev.off()
  
  pdf(file.path(figdir, "L1KCPAll_n=100_countGt100.pdf"), width=8, height=6)
  ggplot() + theme_minimal() +
    geom_line(data=xdf[xdf$compound %in% xdf$compound[xdf$metasize== 50],], aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoBase$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgcpBase$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound metacorrelation, unnormalized")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom", text=element_text(size=16))
  dev.off()
  
  
  #### Example fig: Figure 1B ####
  #xdftop <- xdf[xdf$metasize == 50,]
  xdftop <- xdf[xdf$metasize == 100,]
  xdftop <- xdftop[order(xdftop$meanSim), ]
  
  #examplecps <- xdftop$compound[c(103, 1027, 1952)]
  examplecps <- xdftop$compound[c(37, 378, 713)]
  
  exampledf <- rbind(xdf[xdf$compound %in% examplecps,], 
                     cbind(ndmsoBase$AllDS$cordf, compound="DMSO"), 
                     cbind(nbgcpBase$allCpds$cordf, compound="random"))
  
  pdf(file.path(figdir, "L1KCPExamples_n=100.pdf"), width=8, height=6)
  ggplot(data=exampledf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
    theme_minimal() + geom_ribbon(aes(fill=compound), alpha=0.3) + geom_line(aes(color=compound)) + 
    xlab("Metasignature size") + ylab("Mean Pearson") + coord_cartesian(xlim=c(0, 100)) + 
    ggtitle(sprintf("L1K compound metacorrelation, pan-cell line, unnormalized")) + 
    scale_fill_discrete(breaks=c(examplecps[3:1], "DMSO", "random")) + 
    scale_color_discrete(breaks=c(examplecps[3:1], "DMSO", "random")) + 
    theme(legend.position="bottom", text=element_text(size=16))
  dev.off()
  
  
  #### Distribution figure ####
  refDists <- readRDS(file.path(datadir, "../l1kMeta/simDists/L1KreferenceSimDists.rds"))
  
  # Density plots
  bgsims <- reshape2::melt(refDists$l1kNullCor$allCPDs$simmat)
  
  cpsims <- c()
  metasizes <- rownames(refDists$l1kCPCor[[1]]$simmat)
  for (ii in seq_along(metasizes)){
    cpsims <- rbind(cpsims, data.frame(metasize=metasizes[ii], value=as.numeric(sapply(refDists$l1kCPCor, FUN=function(y) y$simmat[ii,]))))
  }
  cpsims$metasize <- factor(cpsims$metasize, levels = sort(as.numeric(metasizes)))
  
  g1 <- ggplot(bgsims, aes(x=value, group=as.factor(Var1), y=as.factor(Var1), fill=as.factor(Var1))) + 
    geom_density_ridges2(alpha=0.4, scale=1, bandwidth=0.01) + theme_minimal() + xlim(c(-0.4, 1)) + 
    ylab("Metasignature size") + theme(axis.title.x=element_blank()) + # xlab("Pearson's correlation") + 
    guides(fill = guide_legend(title="Metasignature size")) + 
    theme(legend.position="bottom", text=element_text(size=14)) + 
    ggtitle("Metasignatures of Random Compounds") + 
    geom_vline(xintercept = 0, lty=2) + geom_vline(xintercept=1)
  
  g2<- ggplot(cpsims, aes(x=value, group=metasize, y=metasize, fill=metasize)) + 
    geom_density_ridges2(scale=1, bandwidth=0.01) + theme_minimal() + xlim(c(-0.4, 1)) + 
    xlab("Pearson's correlation") + ylab("Metasignature size") + 
    theme(legend.position="none", text=element_text(size=14)) + 
    ggtitle("Compound-specific Metasignatures") + 
    geom_vline(xintercept = 0, lty=2) + geom_vline(xintercept=1)
  
  pdf(file.path(figdir, "L1KPearsonDistExamples.pdf"), width=8, height=6)
  gridExtra::grid.arrange(g1, g2, ncol=1, heights=c(4, 3))
  dev.off()
  
  
  #### Frequency of pos/neg for each gene ####
  posFrac <- readRDS(file.path(datadir, "../l1kMeta/simDists/L1KmetaPosFracs.rds"))
  
  posFracDf <- cbind(reshape2::melt(posFrac$bgMetaPosFrac), reshape2::melt(posFrac$cpMetaPosFrac))
  posFracDf <- posFracDf[, c(4,1,3)]
  colnames(posFracDf) <- c("metasize", "random", "compound")
  posFracDf$geneSymbol <- rep(names(posFrac$bgMetaPosFrac[[1]]), length(posFrac$bgMetaPosFrac))
  
  posFracDf$metasize <- as.numeric(sapply(posFracDf$metasize, FUN=function(x) substr(x, 2, nchar(x))))
  
  pdf(file.path(figdir, "L1KMetaPosFrac.pdf"), width=8, height=6)
  ggplot(posFracDf[posFracDf$metasize %in% c(1,10, 100),], 
         aes(x=random, y=compound, color=as.factor(metasize))) + geom_point(alpha=0.5) + theme_minimal() +
    xlab("Random metasignature positive fraction") + ylab("Compound-specific metasignature positive fraction") +
    theme(legend.position="bottom", text=element_text(size=14)) + guides(color=guide_legend(title="Metasize")) + 
    ggtitle("Fraction of positive landmark genes in metasignatures")
  dev.off()
  
  
  #### Renormalized figures ####
  
  #ndmsoBase <- readRDS(file.path(datadir, "base/dmso", "allDSDMSOMetaSimpearson500x100.rds"))
  #nbgcpBase <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimPearson500x100.rds"))
  
  
  xcp <- readRDS(file.path(datadir, "renormCP/L1KCP", "AllDSMetaSimpearson500x100_compound.rds"))
  xdfcp <- dplyr::bind_rows(lapply(seq_len(length(xcp)), FUN=function(y) cbind(xcp[[y]]$cordf, compound=names(xcp)[y])))
  ndmsoCP <- readRDS(file.path(datadir, "renormCP/dmso/allDSDMSOMetaSimpearson500x100_compound.rds"))
  nbgCP <- readRDS(file.path(datadir, "renormCP/L1Kbg/AllDSBGMetaSimpearson500x100_compound.rds"))
  
  pdf(file.path(figdir, "L1K_RenormCP_CPAll_n=100_countGt100.pdf"), width=8, height=6)
  ggplot() + theme_minimal() +
    geom_line(data=xdfcp[xdfcp$compound %in% xdfcp$compound[xdfcp$metasize== 50],], 
              aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoCP$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgCP$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound metacorrelation: compound normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom", text=element_text(size=14))
  dev.off()
  
  # renorm DMSO:
  xdmso <- readRDS(file.path(datadir, "renormDmso/L1KCP", "AllDSMetaSimpearson500x100_dmso.rds"))
  xdfdmso <- dplyr::bind_rows(lapply(seq_len(length(xdmso)), FUN=function(y) cbind(xdmso[[y]]$cordf, compound=names(xdmso)[y])))
  
  ndmsoDMSO <- readRDS(file.path(datadir, "renormDmso/dmso/allDSDMSOMetaSimpearson500x100_dmso.rds"))
  nbgDMSO <- readRDS(file.path(datadir, "renormDmso/L1Kbg/AllDSBGMetaSimpearson500x100_dmso.rds"))
 
  pdf(file.path(figdir, "L1K_RenormDMSO_CPAll_n=100_countGT100.pdf"), width=8, height=6)
  ggplot() + theme_minimal() +
    geom_line(data=xdfdmso[xdfdmso$compound %in% xdfdmso$compound[xdfdmso$metasize== 50],], 
              aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoDMSO$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgDMSO$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound metacorrelation: DMSO normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom", text=element_text(size=14))
  dev.off()
  
  
  # Read xdf from the Figure 1c section above
  # Plot distribution of similarities for a particular metasignature size for each normalization scheme. 
  
  mymetasize <- 20
  
  xdistDf <- rbind(cbind(xdf[xdf$metasize == mymetasize,], normalize="None"),
                   cbind(xdfcp[xdfcp$metasize == mymetasize,], normalize="Compound"),
                   cbind(xdfdmso[xdfdmso$metasize == mymetasize,], normalize="DMSO"))
  
  pdf(file.path(figdir, sprintf("L1K_RenormSimDists_size=%d.pdf", mymetasize)), width=8, height=6)
  ggplot(xdistDf, aes(x=meanSim, fill=normalize)) + geom_density(alpha=0.3, bw=0.02) + 
    theme_minimal() + xlim(c(0,1)) + xlab("Mean Metasignature Pearson's Correlation") + 
    ylab("Density, bandwidth = 0.02") + 
    ggtitle(sprintf("Metacorrelations for L1000 Compounds, Metasize = %d", mymetasize)) + 
    guides(fill = guide_legend(title="Renormalization")) + 
    theme(legend.position="bottom", text=element_text(size=16))
  dev.off()
  
  
  
  #### Compound renormalization cross validation ####
  
  xvalF <- list.files(file.path(datadir, "renormXval/L1Kbg"))
  
  xvalList <- c()
  for (myxval in xvalF){
    cellid <- strsplit(myxval, "BG")[[1]][1]
    
    xvalds <- readRDS(file.path(datadir, "renormXval/L1Kbg", myxval))
    
    xvalds$allCpds$cordf <- cbind(xvalds$allCpds$cordf, cellid=cellid)
    
    xvalList[cellid] <- xvalds
  }
  
  xvaldf <- dplyr::bind_rows(lapply(xvalList, FUN=function(x) x$cordf))
  
  pdf(file.path(figdir, "L1K_CPRenormXval_metacorr.pdf"), width=8, height=6)
  ggplot(xvaldf, aes(x=metasize, y=meanSim, ymin=meanSim-sdSim, ymax=meanSim+sdSim, color=cellid)) + 
    geom_point(alpha=0.2) + geom_smooth(method="loess") +
    theme_minimal() + scale_x_continuous(trans="log10") + theme(legend.position="bottom") + 
    ylim(c(-0.1, 0.8)) + guides(color = guide_legend(title="Cell line")) + xlab("Metasignature size") + 
    ylab("Mean metasignature Pearson") + theme(text=element_text(size=14)) + 
    ggtitle("Compound renormalization metasignature correlation cross validation")
  dev.off()
  
  #### Power of metasignatures vs random compound metasignatures at alpha = 0.01: ####
  
  # Load data as above:
  xbase <- readRDS(file.path(datadir, "base/L1KCP", "AllDSMetaSimpearson500x100.rds"))
  xdf <- dplyr::bind_rows(lapply(seq_len(length(xbase)), FUN=function(y) cbind(xbase[[y]]$cordf, compound=names(xbase)[y])))
  nbgcpBaseM <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimpearson200x1000.rds"))
 
  xcp <- readRDS(file.path(datadir, "renormCP/L1KCP", "AllDSMetaSimpearson500x100_compound.rds"))
  xdfcp <- dplyr::bind_rows(lapply(seq_len(length(xcp)), FUN=function(y) cbind(xcp[[y]]$cordf, compound=names(xcp)[y])))
  nbgCPM <- readRDS(file.path(datadir, "renormCP/L1Kbg/AllDSBGMetaSimpearson200x1000_compound.rds"))
  
  xdmso <- readRDS(file.path(datadir, "renormDmso/L1KCP", "AllDSMetaSimpearson500x100_dmso.rds"))
  xdfdmso <- dplyr::bind_rows(lapply(seq_len(length(xdmso)), FUN=function(y) cbind(xdmso[[y]]$cordf, compound=names(xdmso)[y])))
  nbgDMSOM <- readRDS(file.path(datadir, "renormDmso/L1Kbg/AllDSBGMetaSimpearson200x1000_dmso.rds"))

  # I have no idea why exactly one element of nbgcpBaseM is NA
  bgBaseThresh <- apply(nbgcpBaseM$allCpds$simmat, 1, FUN=function(x) quantile(x, 0.99, na.rm=1))
  bgCPThresh <- apply(nbgCPM$allCpds$simmat, 1, FUN=function(x) quantile(x, 0.99))
  bgDMSOThresh <- apply(nbgDMSOM$allCpds$simmat, 1, FUN=function(x) quantile(x, 0.99))
  
  
  metaPowerBase <- lapply(seq_along(xbase), FUN=function(x) rowMeans(xbase[[x]]$simmat > bgBaseThresh[seq_len(dim(xbase[[x]]$simmat)[1])]))
  metaPowerCP <- lapply(seq_along(xcp), FUN=function(x) rowMeans(xcp[[x]]$simmat > bgCPThresh[seq_len(dim(xcp[[x]]$simmat)[1])]))
  metaPowerDMSO <- lapply(seq_along(xdmso), FUN=function(x) rowMeans(xdmso[[x]]$simmat > bgDMSOThresh[seq_len(dim(xdmso[[x]]$simmat)[1])]))
  
  ix <- which(sapply(metaPowerBase, length) >= 100)
  powerCurveBase <- sapply(seq_len(200), FUN=function(y) mean(sapply(metaPowerBase[ix], FUN=function(x) x[y]), na.rm=TRUE))
  powerCurveCP <- sapply(seq_len(200), FUN=function(y) mean(sapply(metaPowerCP[ix], FUN=function(x) x[y]), na.rm=TRUE))
  powerCurveDMSO <- sapply(seq_len(200), FUN=function(y) mean(sapply(metaPowerDMSO[ix], FUN=function(x) x[y]), na.rm=TRUE))
  
  saveRDS(list(powerCurveBase=powerCurveBase, 
          powerCurveCP=powerCurveCP,
          powerCurveDMSO=powerCurveDMSO), file=file.path(datadir, "../powerSimulation/powerCurveSummary.rds"))
  
  powerCurvedf <- rbind(data.frame(metasize=seq_along(powerCurveBase), 
                                   metaPower=powerCurveBase,
                                   norm="None"),
                        data.frame(metasize=seq_along(powerCurveCP), 
                                   metaPower=powerCurveCP, 
                                   norm="Compound"), 
                        data.frame(metasize=seq_along(powerCurveDMSO), 
                                   metaPower=powerCurveDMSO,
                                   norm="DMSO"))
  
  pdf(file.path(figdir, "L1K_RenormPower_size=100.pdf"), width=8, height=6)
  ggplot(powerCurvedf, aes(x=metasize, y=metaPower, color=norm)) + geom_line(lwd=2) + theme_minimal() + 
    xlab("Metasignature size") + ylab("Power at alpha = 0.01") + coord_cartesian(xlim=c(0,100)) + 
    geom_hline(yintercept=0.01, lty=2) + guides(color = guide_legend(title="Renormalization")) + 
    theme(legend.position="bottom", text=element_text(size=14)) + 
    ggtitle("Power relative to random compound metasignatures, alpha = 0.01")
  dev.off()
  
  #### Negative control summaries ####
  dmsoSummaryDf <- rbind(cbind(ndmsoBase$AllDS$cordf, normalize="None"), 
                         cbind(ndmsoCP$AllDS$cordf, normalize="Compound"), 
                         cbind(ndmsoDMSO$AllDS$cordf, normalize="DMSO"))
  
  pdf(file.path(figdir, "L1K_DMSORenorm_MetacorCurves.pdf"), width=8, height=6)
  ggplot(dmsoSummaryDf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim + sdSim, color=normalize,
                            fill=normalize)) + geom_line() + geom_ribbon(alpha=0.2) + theme_minimal() + 
    scale_x_continuous(trans="log10") + theme(legend.position="bottom") + xlab("Metasignature size") + 
    ylab("Metasignature Pearson's Correlation") + ggtitle("Distribution of metasignature correlations for DMSO negative controls") + 
    guides(fill=guide_legend(title="Normalization"), color=guide_legend(title="Normalization")) +
    theme(text=element_text(size=14))
  dev.off()
  
  cpbgSummaryDf <- rbind(cbind(nbgcpBase$allCpds$cordf, normalize="None"), 
                         cbind(nbgCP$allCpds$cordf, normalize="Compound"),
                         cbind(nbgDMSO$allCpds$cordf, normalize="DMSO"))
  
  pdf(file.path(figdir, "L1K_BGCPRenorm_MetacorCurves.pdf"), width=8, height=6)
  ggplot(cpbgSummaryDf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim + sdSim, color=normalize,
                            fill=normalize)) + geom_line() + geom_ribbon(alpha=0.2) + theme_minimal() + 
    scale_x_continuous(trans="log10") + theme(legend.position="bottom") + xlab("Metasignature size") + 
    ylab("Metasignature Pearson's Correlation") + ggtitle("Distribution of metasignature correlations for random compounds") + 
    guides(fill=guide_legend(title="Normalization"), color=guide_legend(title="Normalization")) +
    theme(text=element_text(size=14))
  dev.off()
  
  
  #### Offset Power Figures ####
   
  opfiles <- list.files(file.path(datadir, "offsetPower"))
  
  opres <- lapply(opfiles, FUN=function(x) readRDS(file.path(datadir, "offsetPower", x)))
  
  opdf <- dplyr::bind_rows(lapply(opres, FUN=function(x) x$resdf))
  
  pdf(file.path(figdir, "OffsetPowerpLT2_bycell.pdf"), width=8, height=6)
  ggplot(opdf, aes(x=delta, y=pLT01, color=cellid)) + geom_point(alpha=0.5) + geom_smooth(method="loess") + 
    theme_minimal() + xlab("Offset (standard deviations)") + ylab("Power (alpha = 0.01)") + 
    ggtitle("Simulated Offset power: P < 0.01") + ylim(c(0, 0.25)) + theme(text=element_text(size=14))
  dev.off()
  
  pdf(file.path(figdir, "OffsetPowerfdr10_bycell.pdf"), width=8, height=6)
  ggplot(opdf, aes(x=delta, y=fdrLT10, color=cellid)) + geom_point(alpha=0.5) + geom_smooth(method="loess") + 
    theme_minimal() + xlab("Offset (standard deviations)") + ylab("Recall (FDR < 0.10)") + 
    ggtitle("Simulated Offset Recall: FDR < 0.10") + ylim(c(0, 0.4)) + theme(text=element_text(size=14))
  dev.off()
  
  pdf(file.path(figdir, "OffsetPowerfdr25_bycell.pdf"), width=8, height=6)
  ggplot(opdf, aes(x=delta, y=fdrLT25, color=cellid)) + geom_point(alpha=0.5) + geom_smooth(method="loess") + 
    theme_minimal() + xlab("Offset (standard deviations)") + ylab("Recall (FDR < 0.25)") + 
    ggtitle("Simulated Offset Recall: FDR < 0.25") + ylim(c(0, 0.75)) + theme(text=element_text(size=14))
  dev.off()
  
  pdf(file.path(figdir, "OffsetPowerSNR2_bycell.pdf"), width=8, height=6)
  ggplot(opdf, aes(x=delta, y=SNR2, color=cellid)) + geom_point(alpha=0.5) + geom_smooth(method="loess") + 
    theme_minimal() + xlab("Offset (standard deviations)") + ylab("Fraction with SNR > 2") + 
    ggtitle("Simulated Offset Recall: Z-score > 2") + ylim(c(0, 0.35)) + theme(text=element_text(size=14))
  dev.off()
  
}