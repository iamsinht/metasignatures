library(ggplot2)

datadir <- "~/Work/bhk/analysis/metasig/h4h"
figdir <- "~/Work/bhk/analysis/metasig/figures"


#### Unnormalized
# DMSO unnorm
x <- readRDS(file.path(datadir, "dmso", "allCellsDMSOMetaSimpearson250x100.rds"))

dmso0 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

pdf(file.path(figdir, "L1KDMSOByCell_n=100.pdf"), width=10, height=8)
ggplot(dmso0, aes(x=metasize, y=meanSim, color=cell))+ geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation by cell line, N=100 iterations")
dev.off()


# Null: All CP unnorm
BGds0 <- c()
for (myf in list.files(file.path(datadir, "L1KBG"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  
  x <- readRDS(file.path(datadir, "L1KBG", myf))
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


# Compounds: unnorm
pdf(file.path(figdir, "L1KCPByCell_n=100.pdf"), width=10, height=8)

for (myf in list.files(file.path(datadir, "L1KCP"))){
  mycell <- strsplit(myf, split="Meta")[[1]][1]
  print(mycell)
  
  x <- readRDS(file.path(datadir, "L1KCP", myf))
  xdf <- dplyr::bind_rows(lapply(seq_len(length(x)), FUN=function(y) cbind(x[[y]]$cordf, compound=names(x)[y])))
  
  print(ggplot(xdf, aes(x=metasize, y=meanSim, fill=compound)) + geom_line(color="blue") + theme_minimal() + 
    xlab("Metasignature size") + ylab("Mean Pearson") + ggtitle(sprintf("L1K compound correlation, %s, N = 100 iters", mycell)) + 
    scale_x_continuous(trans="log10"))
}
dev.off()


#### DMSO Renorm
# DMSO dmso norm
x <- readRDS(file.path(datadir, "renormDmso/dmso", "allCellsDMSOMetaSimpearson500x100_dmso.rds"))

dmso1 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

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



#### Compound Renorm
# DMSO compound norm
x <- readRDS(file.path(datadir, "renormCP/dmso", "allCellsDMSOMetaSimpearson500x100_compound.rds"))

dmso2 <- dplyr::bind_rows(lapply(seq_along(x), FUN=function(y) cbind(x[[y]]$cordf, cell=names(x)[y])))

pdf(file.path(figdir, "L1KDMSOByCell_n=100_CPNorm.pdf"), width=10, height=8)
ggplot(dmso2, aes(x=metasize, y=meanSim, color=cell))+ geom_line(linewidth=1.5) + theme_minimal() + xlab("Metasignature size") + 
  ylab("Mean Pearson Correlation") + ggtitle("L1K DMSO metasignature correlation by cell line, N=100 iterations")
dev.off()

  
  
  