library(ggplot2)

datadir <- "~/Work/bhk/analysis/metasig/h4h"
figdir <- "~/Work/bhk/analysis/metasig/figures"


#### Unnormalized ####
# DMSO unnorm
x <- readRDS(file.path(datadir, "base/dmso", "allCellsDMSOMetaSimpearson500x100.rds"))

dmso0 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

pdf(file.path(figdir, "L1KDMSOByCell_n=100_x=500.pdf"), width=10, height=8)
ggplot(dmso0, aes(x=metasize, y=meanSim, color=cell))+ geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation by cell line, N=100 iterations")
dev.off()

xall <- readRDS(file.path(datadir, "base/dmso", "allDSDMSOMetaSimpearson500x100.rds"))

ggplot(xall$AllDS$cordf, aes(x=metasize, y=meanSim, ymin=meanSim-sdSim, ymax=meanSim+sdSim)) + 
  geom_ribbon(alpha=0.3) + geom_line() + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation, all cell lines, N = 100 iterations") +
  ylim(c(-0.2, 1))


# Null: All CP unnorm
BGds0 <- c()
for (myf in list.files(file.path(datadir, "base/L1Kbg"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  
  x <- readRDS(file.path(datadir, "base/L1Kbg", myf))
  BGds0[[mycell]] <- cbind(x$allCpds$cordf, cell=mycell)
}

pdf(file.path(figdir, "L1KNullBgByCell_n=100_ribbon.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds0), aes(x=metasize, y=meanSim, fill=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_ribbon(alpha=0.4) + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()

pdf(file.path(figdir, "L1KNullBgByCell_n=100_line.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds0), aes(x=metasize, y=meanSim, color=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()

bgall <- readRDS(file.path(datadir, "base/L1Kbg", "AllDSBGMetaSimPearson500x100.rds"))

ggplot(bgall$allCpds$cordf, aes(x=metasize, y=meanSim, ymin=meanSim-sdSim, ymax=meanSim+sdSim)) + 
  geom_ribbon(alpha=0.3) + geom_line() + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K Background compounds metasignature correlation, all cell lines, N = 100 iterations") +
  ylim(c(-0.2, 1))

# Compounds: unnorm
pdf(file.path(figdir, "L1KCPByCell_n=100.pdf"), width=10, height=8)

for (myf in list.files(file.path(datadir, "base/L1KCP"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  print(mycell)
  
  x <- readRDS(file.path(datadir, "base/L1KCP", myf))
  xdf <- dplyr::bind_rows(lapply(seq_len(length(x)), FUN=function(y) cbind(x[[y]]$cordf, compound=names(x)[y])))
  
  print(ggplot(xdf, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="blue") + theme_minimal() + 
    xlab("Metasignature size") + ylab("Mean Pearson") + ggtitle(sprintf("L1K compound correlation, %s, N = 100 iters", mycell)) + 
    scale_x_continuous(trans="log10"))
}
dev.off()


#### DMSO Renorm ####
# DMSO dmso norm
x <- readRDS(file.path(datadir, "renormDmso/dmso", "allCellsDMSOMetaSimpearson500x100_dmso.rds"))

dmso1 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

# What the actual fuck:
pdf(file.path(figdir, "L1KDMSOByCell_n=100_DMSONorm.pdf"), width=10, height=8)
ggplot(dmso1, aes(x=metasize, y=meanSim, color=cell))+ geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation by cell line, N=100 iterations")
dev.off()


# Null all CP, DMSO norm
BGds1 <- c()
for (myf in list.files(file.path(datadir, "renormDmso/L1KBG"))){
  mycell <- strsplit(myf, split="BGMeta")[[1]][1]
  
  x <- readRDS(file.path(datadir, "renormDmso/L1KBG", myf))
  BGds1[[mycell]] <- cbind(x$allCpds$cordf, cell=mycell)
}

pdf(file.path(figdir, "L1KNullBgByCell_n=100_ribbon_DMSONorm.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds1), aes(x=metasize, y=meanSim, fill=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_ribbon(alpha=0.3) + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()

pdf(file.path(figdir, "L1KNullBgByCell_n=100_line_DMSONorm.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds1), aes(x=metasize, y=meanSim, color=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()


# Compounds, DMSO norm:
pdf(file.path(figdir, "L1KCPByCell_n=100_DMSONorm.pdf"), width=10, height=8)

for (myf in list.files(file.path(datadir, "renormDmso/L1KCP"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  print(mycell)
  
  x <- readRDS(file.path(datadir, "renormDmso/L1KCP", myf))
  xdfdm <- dplyr::bind_rows(lapply(seq_len(length(x)), FUN=function(y) cbind(x[[y]]$cordf, compound=names(x)[y])))
  
  print(ggplot(xdfdm, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="blue") + theme_minimal() + 
          xlab("Metasignature size") + ylab("Mean Pearson") + 
          ggtitle(sprintf("L1K compound correlation, %s, N = 100 iters; DMSO Renormalization", mycell)) + 
          scale_x_continuous(trans="log10"))
}
dev.off()



#### Compound Renorm ####
# DMSO compound norm
x <- readRDS(file.path(datadir, "renormCP/dmso", "allCellsDMSOMetaSimpearson500x100_compound.rds"))

dmso2 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

pdf(file.path(figdir, "L1KDMSOByCell_n=100_CPNorm.pdf"), width=10, height=8)
ggplot(dmso2, aes(x=metasize, y=meanSim, color=cell))+ geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation by cell line, N=100 iterations")
dev.off()


# Null all CP, compound norm:
BGds2 <- c()
for (myf in list.files(file.path(datadir, "renormCP/L1Kbg"))){
  mycell <- strsplit(myf, split="BGMeta")[[1]][1]
  
  x <- readRDS(file.path(datadir, "renormCP/L1Kbg", myf))
  BGds2[[mycell]] <- cbind(x$allCpds$cordf, cell=mycell)
}

pdf(file.path(figdir, "L1KNullBgByCell_n=100_ribbon_CPNorm.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds2), aes(x=metasize, y=meanSim, fill=cell, color=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_ribbon(alpha=0.2) + geom_line() + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()

pdf(file.path(figdir, "L1KNullBgByCell_n=100_line_CPNorm.pdf"), width=10, height=8)
ggplot(dplyr::bind_rows(BGds2), aes(x=metasize, y=meanSim, color=cell, ymin=meanSim - sdSim, ymax=meanSim+sdSim)) + 
  geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + ylab("Mean Pearson") + 
  ggtitle("L1K Background compounds metasignature correlation by cell line, N=100 iters")
dev.off()

  
# Compounds, compound norm:
pdf(file.path(figdir, "L1KCPByCell_n=100_CPNorm.pdf"), width=10, height=8)

for (myf in list.files(file.path(datadir, "renormCP/L1KCP"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  print(mycell)
  
  x <- readRDS(file.path(datadir, "renormCP/L1KCP", myf))
  xdfcp <- dplyr::bind_rows(lapply(seq_len(length(x)), FUN=function(y) cbind(x[[y]]$cordf, compound=names(x)[y])))
  
  print(ggplot(xdfcp, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="blue") + theme_minimal() + 
          xlab("Metasignature size") + ylab("Mean Pearson") + 
          ggtitle(sprintf("L1K compound correlation, %s, N = 100 iters; Compound Renormalization", mycell)) + 
          scale_x_continuous(trans="log10"))
}
dev.off()
 


#### All Compounds, pan-dataset ####
# Base:
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
  ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; base normalization")) +
  scale_x_continuous(trans="log10") + labs(color="Legend") + 
  scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset")
dev.off()

pdf(file.path(figdir, "L1KCPAll_n=100_countGt100.pdf"), width=8, height=6)
ggplot() + theme_minimal() +
  geom_line(data=xdf[xdf$compound %in% xdf$compound[xdf$metasize== 50],], aes(x=metasize, y=meanSim, group=compound, color="Compounds"), alpha=0.3) + 
  geom_line(data=ndmsoBase$AllDS$cordf, aes(x=metasize, y=meanSim, color="DMSO"), linewidth=1) + 
  geom_line(data=nbgcpBase$allCpds$cordf, aes(x=metasize, y=meanSim, color="Random compounds"), linewidth=1) + 
  xlab("Metasignature size") + ylab("Mean Pearson") +
  ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; base normalization")) +
  scale_x_continuous(trans="log10") + labs(color="Legend") + 
  scale_color_manual(values=c("Compounds"="blue", "DMSO"="red", "Random compounds"="black"), name="Dataset")
dev.off()


# renorm CP:
xcp <- readRDS(file.path(datadir, "renormCP/L1KCP", "AllDSMetaSimpearson500x100_compound.rds"))
xdfcp <- dplyr::bind_rows(lapply(seq_len(length(xcp)), FUN=function(y) cbind(xcp[[y]]$cordf, compound=names(xcp)[y])))

pdf(file.path(figdir, "L1KCPAll_n=100_CPNorm.pdf"))
ggplot(xdfcp, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="forestgreen", alpha=0.3) + theme_minimal() +
  xlab("Metasignature size") + ylab("Mean Pearson") +
  ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; Compound normalization")) +
  scale_x_continuous(trans="log10")
dev.off()
  

# renorm DMSO:
xdmso <- readRDS(file.path(datadir, "renormDmso/L1KCP", "AllDSMetaSimpearson500x100_dmso.rds"))
xdfdmso <- dplyr::bind_rows(lapply(seq_len(length(xdmso)), FUN=function(y) cbind(xdmso[[y]]$cordf, compound=names(xdmso)[y])))

pdf(file.path(figdir, "L1KCPAll_n=100_DMSONorm.pdf"))
ggplot(xdfdmso, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="red", alpha=0.3) + theme_minimal() +
  xlab("Metasignature size") + ylab("Mean Pearson") +
  ggtitle(sprintf("L1K compound correlation, pan-cell line, N = 100 iters; DMSO normalization")) +
  scale_x_continuous(trans="log10")
dev.off()



