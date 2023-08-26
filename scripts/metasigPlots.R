library(ggplot2)
library(ggridges)
library(grid)

plotMetasigs <- function(figdir = "~/Work/bhk/analysis/metasig/figures", datadir = "~/Work/bhk/analysis/metasig/h4h"){
  
  #### Figure 1C: ####
  xbase <- readRDS(file.path(datadir, "base/L1KCP", "AllDSMetaSimpearson500x100.rds"))
  xdf <- dplyr::bind_rows(lapply(seq_len(length(xbase)), FUN=function(y) cbind(xbase[[y]]$cordf, compound=names(xbase)[y])))
  
  ndmsoBase <- readRDS(file.path(datadir, "base/dmso", "allDSDMSOMetaSimpearson500x100.rds"))
  nbgcpBase <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimPearson500x100.rds"))
  
  pdf(file.path(figdir, "L1KCPAll_n=100.pdf"), width=8, height=7)
  ggplot() + theme_minimal() +
    geom_line(data=xdf, aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoBase$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgcpBase$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; base normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom")
  dev.off()
  
  pdf(file.path(figdir, "L1KCPAll_n=100_countGt100.pdf"), width=8, height=7)
  ggplot() + theme_minimal() +
    geom_line(data=xdf[xdf$compound %in% xdf$compound[xdf$metasize== 50],], aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoBase$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgcpBase$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; base normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom")
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
  
  pdf(file.path(figdir, "L1KCPExamples_n=100.pdf"), width=8, height=7)
  ggplot(data=exampledf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
    theme_minimal() + geom_ribbon(aes(fill=compound), alpha=0.3) + geom_line(aes(color=compound)) + 
    xlab("Metasignature size") + ylab("Mean Pearson") + coord_cartesian(xlim=c(0, 100)) + 
    ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; base normalization")) + 
    scale_fill_discrete(breaks=c(examplecps[3:1], "DMSO", "random")) + 
    scale_color_discrete(breaks=c(examplecps[3:1], "DMSO", "random")) + theme(legend.position="bottom")
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
    xlab("Pearson's correlation") + ylab("Metasignature size") + 
    guides(fill = guide_legend(title="Metasignature size")) + theme(legend.position="bottom") + 
    ggtitle("Metasignatures of Random Compounds") + 
    geom_vline(xintercept = 0, lty=2) + geom_vline(xintercept=1)
  
  g2<- ggplot(cpsims, aes(x=value, group=metasize, y=metasize, fill=metasize)) + 
    geom_density_ridges2(scale=1, bandwidth=0.01) + theme_minimal() + xlim(c(-0.4, 1)) + 
    xlab("Pearson's correlation") + ylab("Metasignature size") + theme(legend.position="none") + 
    ggtitle("Compound-specific Metasignatures") + 
    geom_vline(xintercept = 0, lty=2) + geom_vline(xintercept=1)
  
  pdf(file.path(figdir, "L1KPearsonDistExamples.pdf"), width=8, height=7)
  gridExtra::grid.arrange(g1, g2, ncol=1, heights=c(4, 3))
  dev.off()
  
  
  #### Frequency of pos/neg for each gene ####
  posFrac <- readRDS(file.path(datadir, "../l1kMeta/simDists/L1KmetaPosFracs.rds"))
  
  posFracDf <- cbind(reshape2::melt(posFrac$bgMetaPosFrac), reshape2::melt(posFrac$cpMetaPosFrac))
  posFracDf <- posFracDf[, c(4,1,3)]
  colnames(posFracDf) <- c("metasize", "random", "compound")
  posFracDf$geneSymbol <- rep(names(posFrac$bgMetaPosFrac[[1]]), length(posFrac$bgMetaPosFrac))
  
  posFracDf$metasize <- as.numeric(sapply(posFracDf$metasize, FUN=function(x) substr(x, 2, nchar(x))))
  
  pdf(file.path(figdir, "L1KMetaPosFrac.pdf"), width=8, height=7)
  ggplot(posFracDf[posFracDf$metasize %in% c(1,10, 100),], 
         aes(x=random, y=compound, color=as.factor(metasize))) + geom_point(alpha=0.5) + theme_minimal() +
    xlab("Random metasignature positive fraction") + ylab("Compound-specific metasignature positive fraction") +
    theme(legend.position="bottom") + guides(color=guide_legend(title="Metasize"))
  dev.off()
  
  
  #### Renormalized figures ####
  
  #ndmsoBase <- readRDS(file.path(datadir, "base/dmso", "allDSDMSOMetaSimpearson500x100.rds"))
  #nbgcpBase <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimPearson500x100.rds"))
  
  
  xcp <- readRDS(file.path(datadir, "renormCP/L1KCP", "AllDSMetaSimpearson500x100_compound.rds"))
  xdfcp <- dplyr::bind_rows(lapply(seq_len(length(xcp)), FUN=function(y) cbind(xcp[[y]]$cordf, compound=names(xcp)[y])))
  ndmsoCP <- readRDS(file.path(datadir, "renormCP/dmso/allDSDMSOMetaSimpearson500x100_compound.rds"))
  nbgCP <- readRDS(file.path(datadir, "renormCP/L1Kbg/AllDSBGMetaSimpearson500x100_compound.rds"))
  
  pdf(file.path(figdir, "L1K_RenormCP_CPAll_n=100_countGt100.pdf"), width=8, height=7)
  ggplot() + theme_minimal() +
    geom_line(data=xdfcp[xdfcp$compound %in% xdfcp$compound[xdfcp$metasize== 50],], 
              aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoCP$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgCP$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; compound normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom")
  dev.off()
  
  # renorm DMSO:
  xdmso <- readRDS(file.path(datadir, "renormDmso/L1KCP", "AllDSMetaSimpearson500x100_dmso.rds"))
  xdfdmso <- dplyr::bind_rows(lapply(seq_len(length(xdmso)), FUN=function(y) cbind(xdmso[[y]]$cordf, compound=names(xdmso)[y])))
  
  ndmsoDMSO <- readRDS(file.path(datadir, "renormDmso/dmso/allDSDMSOMetaSimpearson500x100_dmso.rds"))
  nbgDMSO <- readRDS(file.path(datadir, "renormDmso/L1Kbg/AllDSBGMetaSimpearson500x100_dmso.rds"))
 
  pdf(file.path(figdir, "L1K_RenormDMSO_CPAll_n=100_countGT100.pdf"), width=8, height=7)
  ggplot() + theme_minimal() +
    geom_line(data=xdfdmso[xdfdmso$compound %in% xdfdmso$compound[xdfdmso$metasize== 50],], 
              aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
    geom_line(data=ndmsoDMSO$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
    geom_line(data=nbgDMSO$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
    xlab("Metasignature size") + ylab("Mean Pearson") +
    ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; compound normalization")) +
    scale_x_continuous(trans="log10") + labs(color="Legend") + 
    scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset") +
    theme(legend.position="bottom")
  dev.off()
  
  
  # Read xdf from the Figure 1c section above
  # Plot distribution of similarities for a particular metasignature size for each normalization scheme. 
  
  mymetasize <- 20
  
  xdistDf <- rbind(cbind(xdf[xdf$metasize == mymetasize,], normalize="None"),
                   cbind(xdfcp[xdfcp$metasize == mymetasize,], normalize="Compound"),
                   cbind(xdfdmso[xdfdmso$metasize == mymetasize,], normalize="DMSO"))
  
  pdf(file.path(figdir, sprintf("L1K_RenormSimDists_size=%d.pdf", mymetasize)), width=8, height=7)
  ggplot(xdistDf, aes(x=meanSim, fill=normalize)) + geom_density(alpha=0.3, bw=0.02) + 
    theme_minimal() + xlim(c(0,1)) + xlab("Mean Metasignature Pearson's Correlation") + 
    ylab("Density, bandwidth = 0.02") + 
    ggtitle(sprintf("Distribution of metasignature correlations for L1000 Compounds, Metasize = %d", mymetasize)) + 
    guides(fill = guide_legend(title="Renormalization")) + theme(legend.position="bottom")
  dev.off()
  
  
  
  # Negative control summaries
  dmsoSummaryDf <- rbind(cbind(ndmsoBase$AllDS$cordf, normalize="None"), 
                         cbind(ndmsoCP$AllDS$cordf, normalize="Compound"), 
                         cbind(ndmsoDMSO$AllDS$cordf, normalize="DMSO"))
  
  pdf(file.path(figdir, "L1K_DMSORenorm_MetacorCurves.pdf"), width=8, height=7)
  ggplot(dmsoSummaryDf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim + sdSim, color=normalize,
                            fill=normalize)) + geom_line() + geom_ribbon(alpha=0.2) + theme_minimal() + 
    scale_x_continuous(trans="log10") + theme(legend.position="bottom") + xlab("Metasignature size") + 
    ylab("Metasignature Pearson's Correlation") + ggtitle("Distribution of metasignature correlations for DMSO negative controls") + 
    guides(fill=guide_legend(title="Normalization"), color=guide_legend(title="Normalization")) 
  dev.off()
  
  cpbgSummaryDf <- rbind(cbind(nbgcpBase$allCpds$cordf, normalize="None"), 
                         cbind(nbgCP$allCpds$cordf, normalize="Compound"),
                         cbind(nbgDMSO$allCpds$cordf, normalize="DMSO"))
  
  pdf(file.path(figdir, "L1K_BGCPRenorm_MetacorCurves.pdf"), width=8, height=7)
  ggplot(cpbgSummaryDf, aes(x=metasize, y=meanSim, ymin=meanSim - sdSim, ymax=meanSim + sdSim, color=normalize,
                            fill=normalize)) + geom_line() + geom_ribbon(alpha=0.2) + theme_minimal() + 
    scale_x_continuous(trans="log10") + theme(legend.position="bottom") + xlab("Metasignature size") + 
    ylab("Metasignature Pearson's Correlation") + ggtitle("Distribution of metasignature correlations for random compounds") + 
    guides(fill=guide_legend(title="Normalization"), color=guide_legend(title="Normalization")) 
  dev.off()
   
}