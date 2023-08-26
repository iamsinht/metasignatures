library(metasignatures)
library(cmapR)
library(perturbKit)

# While the metasignatures and perturbKit packages contain the bulk of the functionality,
# the data generation calls are not in general stored, as they are submitted as jobs. 
# This script is designed to contain all the data generation calls, though in practice,
# many of these function calls were run as jobs on an HPC. 

datapath <- "~/Work/bhk/data/l1k/2020"
datadir <- "~/Work/bhk/analysis/metasig/h4h"
figdir <- "~/Work/bhk/analysis/metasig/figures"

# Calculate Pearson correlation distributions (Figure 1d)
x <- l1kSimDists(datapath, datapath, file.path(datadir, "../l1kMeta/simDists"))

# Calculate positive fraction (Figure 1e)
ret <- l1kMetaPosFrac(dspath=datapath, metapath=datapath, outpath=file.path(figdir, "../l1kMeta/simDists/"))


# Calculate moments

# base

# Metasignatures
l1kmeta <- perturbKit::read_l1k_meta(datapath, version=2020)
attach(l1kmeta)
ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)

ix <- sample(seq_along(ds@cid), 100000)
meta10 <- getMetasigs(ds@mat[, ix], k=10, return2=0, returnk=1)
meta100 <- getMetasigs(ds@mat[, ix], k=100, return2=0, returnk=1)

moment10 <- getMoments(meta10, file.path(datadir, "../moments/metasigs/L1K_all_k=10_LM.rds"))
moment100 <- getMoments(meta100, file.path(datadir, "../moments/metasigs/L1K_all_k=100_LM.rds"))



# Power simulations
l1kmeta <- perturbKit::read_l1k_meta(datapath, version=2020)
outdir <- "~/Work/bhk/analysis/metasig/powerSimulation"

topcells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE))[1:30]
offpower <- c()

for (mycell in topcells[10:1]){
  print(mycell)
  offpower[[mycell]] <- offsetPower(cellline=mycell, datapath=datapath, metapath=datapath, 
                        outpath=file.path(outdir, sprintf("%sOffSetPower_iter=10.rds", mycell)), iter=10, l1kmeta=l1kmeta)
}

offpowerdf <- dplyr::bind_rows(lapply(offpower, FUN=function(x) x$resdf))

ggplot(offpowerdf, aes(x=delta, y=fdrLT25, group=cellid, color=cellid)) + geom_point() + geom_smooth(method="loess")