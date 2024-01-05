# library(perturbKit)
# library(cmapR)


#' l1kMetaCalc
#' 
#' Compute correlation of metasignatures for L1000 cell lines
#' @param datapath Path to directory containing level5 ds 
#' @param metapath Path to directory containing metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param renorm Whether to renormalize the data. Default 0. One of c(0, "compound", "dmso")
#' @param metric Metric to use in getMetaSims. Default is Pearson. 
#' @param normpath Path to dataset input for normalization (optional). Default c()
#' @param bycell Logical, whether to run by cell or using the entire dataset. 
#' 
#' @returns Saves output, doesn't return anything. 
#' @export
l1kMetaCalc <- function(datapath=".", metapath=".", outpath=".", 
                        kmax=100, iter=100, renorm="0", metric="pearson",
                        normpath = c(), bycell=TRUE){

  renorm <- match.arg(renorm, c("compound", "dmso", 0))
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  if (bycell){
    # Get the 25 most assayed L1K cell lines
    mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
    
    for (mycell in mycells){
      print(mycell)
      ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                              cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], 
                              rid = landmarks$pr_gene_id)
      
      if (renorm != "0"){
        if (renorm == "compound"){
          ds@mat <- renormData(ds@mat, ds@mat, method="center")
        } else if (renorm == "dmso"){
          normsigs <- siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_iname == "DMSO"]
          normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
          
          ds@mat <- renormData(normds@mat, ds@mat, method="center")
        }
      }
      
      mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
      
      l1kMetaCor <- getMetaSimDs(ds@mat, mysigs$pert_iname, kmax=kmax, iter=iter, metric=metric)
      
      outfile <- sprintf("%sMetaSim%s%dx%d.rds", mycell, metric, kmax, iter)
      if (renorm != "0"){
        outfile <- sprintf("%sMetaSim%s%dx%d_%s.rds", mycell, metric, kmax, iter, renorm)
      }
      
      saveRDS(l1kMetaCor, file.path(outpath, outfile))
    }
  } else {
    
    print("Using entire dataset irrespective of cell lines for metasignature calculation.")
    
    ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                            cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], 
                            rid = landmarks$pr_gene_id)
    
    if (renorm != "0"){
      if (renorm == "compound"){
        ds@mat <- renormData(ds@mat, ds@mat, method="center")
      } else if (renorm == "dmso"){
        normsigs <- siginfo$sig_id[siginfo$pert_iname == "DMSO"]
        normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
        
        ds@mat <- renormData(normds@mat, ds@mat, method="center")
      }
    }
    
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    l1kMetaCor <- getMetaSimDs(ds@mat, mysigs$pert_iname, kmax=kmax, iter=iter, metric=metric)
    
    outfile <- sprintf("AllDSMetaSim%s%dx%d.rds", metric, kmax, iter)
    if (renorm != "0"){
      outfile <- sprintf("AllDSMetaSim%s%dx%d_%s.rds", metric, kmax, iter, renorm)
    }
    
    saveRDS(l1kMetaCor, file.path(outpath, outfile))
  }
}


#' l1kBgCalc
#' 
#' Compute correlation of background metasignatures for L1000 cell lines
#' @param datapath Path to directory containing level5 ds 
#' @param metapath Path to directory containing metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param renorm Whether to renormalize the data. Default 0. One of c(0, "compound", "dmso")
#' @param metric Metric to use in getMetaSims. Default is Pearson. 
#' @param normpath Path to dataset input for normalization (optional). Default c()
#' @param bycell Logical, whether to run by cell or using the entire dataset. 
#' 
#' @returns Saves output, doesn't return anything. 
#' @export
l1kBgCalc <- function(datapath=".", metapath=".", outpath=".", 
                      kmax=100, iter=100, renorm=0, metric="pearson", 
                      normpath=c(), bycell=TRUE){
  
  renorm <- match.arg(renorm, c("compound", "dmso", "xval", 0))
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  if (bycell){
    # Get the 25 most assayed L1K cell lines
    mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
    
    for (mycell in mycells){
      print(mycell)
      ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                              cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], 
                              rid = landmarks$pr_gene_id)
      
      if (renorm != "0"){
        if (renorm == "compound"){
          ds@mat <- renormData(ds@mat, ds@mat, method="center")
        } else if (renorm == "xval"){
          ix <- sample(length(ds@cid), floor(length(ds@cid)/2))

          ds1 <- cmapR::subset_gct(ds, cid=ix)
          ds2 <- cmapR::subset_gct(ds, cid=setdiff(seq(length(ds@cid)), ix))
          
          # This returns ds2 centered using ds1's mean, though the datasets are interchangeable
          ds2@mat <- renormData(ds1@mat, ds2@mat, method="center")
          ds <- ds2
          
        }else if (renorm == "dmso"){
          normsigs <- siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_iname == "DMSO"]
          normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
          
          ds@mat <- renormData(normds@mat, ds@mat, method="center")
        }
      }
      
      mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
      
      l1kMetaCor <- getMetaSimDs(ds@mat, rep("allCpds", dim(ds@mat)[2]),  kmax=kmax, iter=iter, metric=metric)
      
      outfile <- sprintf("%sBGMetaSim%s%dx%d.rds", mycell, metric, kmax, iter)
      if (renorm != "0"){
        outfile <- sprintf("%sBGMetaSim%s%dx%d_%s.rds", mycell, metric, kmax, iter, renorm)
      }
      
      saveRDS(l1kMetaCor, file.path(outpath, outfile))
    }
  } else {
    print("Using entire dataset irrespective of cell lines for metasignature background calculation.")
    
    ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                            cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], 
                            rid = landmarks$pr_gene_id)
    
    if (renorm != "0"){
      if (renorm == "compound"){
        ds@mat <- renormData(ds@mat, ds@mat, method="center")
      } else if (renorm == "xval"){
        ix <- sample(length(ds@cid), floor(length(ds@cid)/2))
        
        ds1 <- cmapR::subset_gct(ds, cid=ix)
        ds2 <- cmapR::subset_gct(ds, cid=setdiff(seq(length(ds@cid)), ix))
        
        # This returns ds2 centered using ds1's mean, though the datasets are interchangeable
        ds@mat <- renormData(ds1@mat, ds2@mat, method="center")
        
      } else if (renorm == "dmso"){
        normsigs <- siginfo$sig_id[siginfo$pert_iname == "DMSO"]
        normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
        
        ds@mat <- renormData(normds@mat, ds@mat, method="center")
      }
    }
    
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    l1kMetaCor <- getMetaSimDs(ds@mat, rep("allCpds", dim(ds@mat)[2]),  kmax=kmax, iter=iter, metric=metric)
    
    outfile <- sprintf("%sBGMetaSim%s%dx%d.rds", "AllDS", metric, kmax, iter)
    if (renorm != "0"){
      outfile <- sprintf("%sBGMetaSim%s%dx%d_%s.rds", "AllDS", metric, kmax, iter, renorm)
    }
    
    saveRDS(l1kMetaCor, file.path(outpath, outfile))
  }
  
}


#' l1kNullCalc
#' 
#' Computes correlation of metasignatures of negative control DMSOs
#' @param dspath Path to level5 control dataset
#' @param metapath Path to directory containing L1000 metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param metric Metric to use in getMetaSims, default is "pearson"
#' @param renorm Whether to renormalize the data. Default 0. Options (0, "compound", "dmso")
#' @param normpath Path to dataset input for compound normalization (optional). Default c()
#' @param bycell Logical, whether to run by cell or using the entire dataset. 
#' 
#' @returns Saves output, returns list of getMetaSimDs objects. 
#' @export
l1kNullCalc <- function(dspath=".", metapath=".", outpath=".", 
                        kmax=100, iter=100, metric="pearson", renorm="0", 
                        normpath=c(), bycell=TRUE){

  renorm <- match.arg(renorm, c("compound", "dmso", "0"))
    
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks

  ds <- cmapR::parse_gctx(dspath, rid=landmarks$pr_gene_id)
  dmsores <- list()
  
  if (bycell == TRUE){  
    # Get the 25 most assayed L1K cell lines
    mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    for (mycell in mycells){
      print(mycell)
      ds1 <- cmapR::subset_gct(ds, cid=mysigs$sig_id[mysigs$cell_id == mycell & mysigs$pert_iname == "DMSO"])
      
      if (renorm != "0"){
        if (renorm == "dmso"){
          ds1@mat <- renormData(ds1@mat, ds1@mat, method="center")
        } else if (renorm == "compound"){
          normsigs <- siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"]
          normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
          
          ds1@mat <- renormData(normds@mat, ds1@mat, method="center")
        }
      }
      
      dmsoMetaCor <- getMetaSimDs(ds1@mat, rep("dmso", dim(ds1@mat)[2]), kmax=kmax, iter=iter, metric=metric)
      
      dmsores[length(dmsores)+1] <- dmsoMetaCor
    }
    
    names(dmsores) <- mycells
    
    outfile <- sprintf("allCellsDMSOMetaSim%s%dx%d.rds", metric, kmax, iter)
    if (renorm != "0"){
      outfile <- sprintf("allCellsDMSOMetaSim%s%dx%d_%s.rds", metric, kmax, iter, renorm)
    }
    
    saveRDS(dmsores, file.path(outpath, outfile))
  
    #return(dmsores)
  } else {
    print("Using entire dataset irrespective of cell lines for DMSO metasignature calculation.")
    
    if (renorm != "0"){
      if (renorm == "dmso"){
        ds@mat <- renormData(ds@mat, ds@mat, method="center")
      } else if (renorm == "compound"){
        normsigs <- siginfo$sig_id[siginfo$pert_type == "trt_cp"]
        normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
        
        ds@mat <- renormData(normds@mat, ds@mat, method="center")
      }
    }
    
    dmsoMetaCor <- getMetaSimDs(ds@mat, rep("dmso", dim(ds@mat)[2]), kmax=kmax, iter=iter, metric=metric)
    
    dmsores[length(dmsores)+1] <- dmsoMetaCor
    names(dmsores) <- c("AllDS")
    
    outfile <- sprintf("allDSDMSOMetaSim%s%dx%d.rds", metric, kmax, iter)
    if (renorm != "0"){
      outfile <- sprintf("allDSDMSOMetaSim%s%dx%d_%s.rds", metric, kmax, iter, renorm)
    }
    
    saveRDS(dmsores, file.path(outpath, outfile))
  }
}


#' l1kSimDists
#' 
#' Compute similarity distributions and fraction of gene values that are positive for L1000
#' data to illustrate the metasignature anomaly.
#' 
#' @inheritParams l1kNullCalc
#' 
#' @returns list of null getMetaSimDS and compound-specific getMetaSimDS outputs.
l1kSimDists <- function(dspath=".", metapath=".", outpath="."){
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                          cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], 
                          rid = landmarks$pr_gene_id)

  pertcount <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
  topperts <- names(pertcount)[pertcount > 200]
  
  topperts <- sample(topperts, 100)
  
  dsTop <- cmapR::subset_gct(ds, cid=siginfo$sig_id[which(siginfo$pert_iname %in% topperts & siginfo$pert_type == "trt_cp")])
  sigsTop <- siginfo[match(dsTop@cid, siginfo$sig_id),]
    
  l1kNullCor <- getMetaSimDs(ds@mat, rep("allCPDs", dim(ds@mat)[2]), kmax=100, iter=10000, seqvals=c(3, 10, 30, 100))
  
  
  l1kCPCor <- getMetaSimDs(dsTop@mat, sigsTop$pert_iname, kmax=100, iter=100, seqvals=c(3, 10, 30, 100))
  
  saveRDS(list(l1kNullCor=l1kNullCor, l1kCPCor=l1kCPCor), file.path(outpath, sprintf("L1KreferenceSimDists.rds")))
  return(list(l1kNullCor=l1kNullCor, l1kCPCor=l1kCPCor))
}


#' l1kMetaSigs
#' 
#' Extract metasignatures of various values to illustrate the tendencies of genes. 
#' 
#' @inheritParams l1kNullCalc
l1kMetaPosFrac <- function(dspath, metapath, outpath="."){
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                          cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], 
                          rid = landmarks$pr_gene_id)
  
  pertcount <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
  topperts <- names(pertcount)[pertcount > 200]
  
  print("Running background")
  ix <- sample(seq_along(ds@cid), 100000)
  bgMetaPosFrac <- list(n1=rowMeans(ds@mat[, ix] > 0),  
                 n3=rowMeans(getMetasigs(ds@mat[, ix], k=3, returnk=1, return2=0) > 0),
                 n10=rowMeans(getMetasigs(ds@mat[, ix], k=10, returnk=1, return2=0) > 0),
                 n30=rowMeans(getMetasigs(ds@mat[, ix], k=30, returnk=1, return2=0) > 0),
                 n100=rowMeans(getMetasigs(ds@mat[, ix], k=100, returnk=1, return2=0) > 0))
  
  print("Running compound-specific")
  print(system.time(cpmetasigs <- sapply(topperts, FUN=function(x) 
    getMetasigs(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname == x]], k=10, return2=0, returnk=1))))
  
  cpMetaPosFrac <- list(n1=rowMeans(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname %in% topperts]] > 0), 
                 n3=rowMeans(Reduce(cbind, sapply(topperts, FUN=function(x) 
                   getMetasigs(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname == x]], k=3, return2=0, returnk=1))) > 0),
                 n10=rowMeans(Reduce(cbind, sapply(topperts, FUN=function(x) 
                   getMetasigs(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname == x]], k=10, return2=0, returnk=1))) > 0),
                 n30=rowMeans(Reduce(cbind, sapply(topperts, FUN=function(x) 
                   getMetasigs(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname == x]], k=30, return2=0, returnk=1))) > 0),
                 n100=rowMeans(Reduce(cbind, sapply(topperts, FUN=function(x) 
                   getMetasigs(ds@mat[, ds@cid %in% siginfo$sig_id[siginfo$pert_iname == x]], k=100, return2=0, returnk=1))) > 0))
  
  ret <- list(bgMetaPosFrac=bgMetaPosFrac, cpMetaPosFrac=cpMetaPosFrac)
  saveRDS(ret, file = file.path(outpath, "L1KmetaPosFracs.rds"))
  
  return(ret)    
}



# The objective is to generate compound-specific metasignatures and compute the similarities
# between signatures of different compounds. The hypothesis is that even with normalization, we
# will still see an increase in cross-compound metasignature correlation with increasing
# metasize. 
l1kMetasigCrossSim <- function(dspath=".", metapath=".", outpath=".", kmax=100, 
                          metric="pearson", renorm="0", normpath=c(), bycell=TRUE){
  
  renorm <- match.arg(renorm, c("compound", "dmso", 0))
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  # By cell is not implemented currently; the code is from l1kMetaCalc
  if (bycell){
    return(0)
    # Get the 25 most assayed L1K cell lines
    mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
    
    for (mycell in mycells){
      print(mycell)
      ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                              cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], 
                              rid = landmarks$pr_gene_id)
      
      if (renorm != "0"){
        if (renorm == "compound"){
          ds@mat <- renormData(ds@mat, ds@mat, method="center")
        } else if (renorm == "dmso"){
          normsigs <- siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_iname == "DMSO"]
          normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
          
          ds@mat <- renormData(normds@mat, ds@mat, method="center")
        }
      }
      
      mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
      
      l1kMetaCor <- getMetaSimDs(ds@mat, mysigs$pert_iname, kmax=kmax, iter=iter, metric=metric)
      
      outfile <- sprintf("%sMetaSim%s%dx%d.rds", mycell, metric, kmax, iter)
      if (renorm != "0"){
        outfile <- sprintf("%sMetaSim%s%dx%d_%s.rds", mycell, metric, kmax, iter, renorm)
      }
      
      saveRDS(l1kMetaCor, file.path(outpath, outfile))
    }
  } else {
    
    print("Computing metasignatures on entire dataset")
    
    ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                            cid=siginfo$sig_id[siginfo$pert_type == "trt_cp"], 
                            rid = landmarks$pr_gene_id)
    
    if (renorm != "0"){
      if (renorm == "compound"){
        ds@mat <- renormData(ds@mat, ds@mat, method="center")
      } else if (renorm == "dmso"){
        normsigs <- siginfo$sig_id[siginfo$pert_iname == "DMSO"]
        normds <- cmapR::parse_gctx(normpath, cid = normsigs, rid = landmarks$pr_gene_id)
        
        ds@mat <- renormData(normds@mat, ds@mat, method="center")
      }
    }
    
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    pertcount <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
    topperts <- names(pertcount)[pertcount > 200]
    
    #grab 250 perts
    set.seed(100)
    topperts <- sample(topperts, 250)

    metaCrossSim <- c()
    kvals <- seq(10, 100, 10)
    
    for (k in kvals){
      print(k)
      pertsigs <- c()

      print(sprintf("Calculating metasignatures, k = %d", k))  
      for (mypert in topperts){
        x <- getMetasigs(ds@mat[, which(mysigs$pert_iname == mypert)], k=k, return2=0, returnk=1)
        pertsigs[[mypert]] <- x[, 1:min(dim(x)[2], 5)]
      }
      
      metablock <- Reduce(cbind, pertsigs)
      
      print("Calculating correlations")
      metacor <- cor(metablock, method="pearson")
      
      pertlengths <- as.numeric(sapply(pertsigs, FUN=function(x) dim(x)[2]))
      pertident <- as.numeric(unlist(sapply(seq_along(pertlengths), FUN=function(x) rep(x, pertlengths[x]))))
      mymask <- outer(pertident, pertident, '==')
    
      metaCrossSim[[length(metaCrossSim)+1]] <- list(metacor=metacor, 
                                                     mymask=mymask,
                                                     metasize=k)
    }
    
    outfile <- sprintf("metaCrossSim_Base.rds")
    if (renorm != "0"){
      outfile <- sprintf("metaCrossSim_%s.rds", renorm)
    }
    
    saveRDS(metaCrossSim, file.path(outpath, outfile))
  }
}
