library(ggplot2)
library(ggridges)
library(grid)

plotMetasigs <- function(figdir = "~/Work/bhk/analysis/metasig/figures"){
  
  # Figure 1C:
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
  
  
  # Example fig: Figure 1B
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
  
  
  # Distribution figure
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
  
  
  ### Frequency of pos/neg for each gene
  
  
}