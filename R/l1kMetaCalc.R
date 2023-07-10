# library(perturbKit)
# library(cmapR)


#datapath <- "~/Work/bhk/data/l1k/2020"
#outpath <- "~/Work/bhk/analysis/metasig/l1kMeta"

#' l1kMetaCalc
#' Compute correlation of metasignatures for L1000 cell lines
#' @param datapath Path to directory containing level5 ds 
#' @param metapath Path to directory containing metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param renorm Whether to renormalize the data, currently not implemented. Default 0. 
#' @param metric Metric to use in getMetaSims. Default is Pearson. 
#' 
#' @returns Saves output, doesn't return anything. 
#' @export
l1kMetaCalc <- function(datapath=".", metapath=".", outpath=".", kmax=100, iter=100, renorm=0, metric="pearson"){

  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  # Get the 25 most assayed L1K cell lines
  mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
  
  for (mycell in mycells){
    print(mycell)
    ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                            cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], 
                            rid = landmarks$pr_gene_id)
    
    if (renorm != 0){
      
    }
    
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    l1kMetaCor <- getMetaSimDs(ds@mat, mysigs$pert_iname, kmax=kmax, iter=iter, metric=metric)
    
    saveRDS(l1kMetaCor, file.path(outpath, sprintf("%sMetaSim%s%dx%d.rds", mycell, metric, kmax, iter)))
  }

}


#' l1kBgCalc
#' Compute correlation of background metasignatures for L1000 cell lines
#' @param datapath Path to directory containing level5 ds 
#' @param metapath Path to directory containing metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param renorm Whether to renormalize the data, currently not implemented. Default 0. 
#' @param metric Metric to use in getMetaSims. Default is Pearson. 
#' 
#' @returns Saves output, doesn't return anything. 
#' @export
l1kBgCalc <- function(datapath=".", metapath=".", outpath=".", kmax=100, iter=100, renorm=0, metric="pearson"){
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks
  
  # Get the 25 most assayed L1K cell lines
  mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
  
  for (mycell in mycells){
    print(mycell)
    ds <- cmapR::parse_gctx(perturbKit::get_level5_ds(datapath, mypattern="trt_cp"), 
                            cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], 
                            rid = landmarks$pr_gene_id)
    
    if (renorm != 0){
      
    }
    
    mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
    
    l1kMetaCor <- getMetaSimDs(ds@mat, rep("allCpds", dim(ds@mat)[2]),  kmax=kmax, iter=iter, metric=metric)
    
    saveRDS(l1kMetaCor, file.path(outpath, sprintf("%sMetaSim%s%dx%d.rds", mycell, metric, kmax, iter)))
  }
  
}


#' l1kNullCalc
#' Computes correlation of metasignatures of negative control DMSOs
#' @param dspath Path to level5 control dataset
#' @param metapath Path to directory containing L1000 metadata
#' @param outpath Path for output files
#' @param kmax Maximum metasignature size to use
#' @param iter Number of iterations for each value of k
#' @param metric Metric to use in getMetaSims, default is "pearson"
#' @param renorm Whether to renormalize the data, currently not implemented. Default 0. 
#' 
#' @returns Saves output, returns list of getMetaSimDs objects. 
#' @export
l1kNullCalc <- function(dspath=".", metapath=".", outpath=".", kmax=100, iter=100, metric="pearson", renorm=0){
  
  l1kmeta <- perturbKit::read_l1k_meta(metapath, version=2020)
  siginfo <- l1kmeta$siginfo
  landmarks <- l1kmeta$landmarks

  ds <- cmapR::parse_gctx(dspath, rid=landmarks$pr_gene_id)
  
  # Get the 25 most assayed L1K cell lines
  mycells <- names(sort(table(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), decreasing=TRUE)[1:25])
  mysigs <- siginfo[match(ds@cid, siginfo$sig_id),]
  
  dmsores <- list()
  
  for (mycell in mycells){
    print(mycell)
    ds1 <- cmapR::subset_gct(ds, cid=mysigs$sig_id[mysigs$cell_id == mycell & mysigs$pert_iname == "DMSO"])
    
    if (renorm != 0){
      
    }
    
    dmsoMetaCor <- getMetaSimDs(ds1@mat, rep("dmso", dim(ds1@mat)[2]), kmax=kmax, iter=iter, metric=metric)
    
    dmsores[length(dmsores)+1] <- dmsoMetaCor
  }
  
  names(dmsores) <- mycells
  
  saveRDS(dmsores, file.path(outpath, sprintf("allCellsDMSOMetaSim%s%dx%d.rds", metric, kmax, iter)))
  detach()
  
  return(dmsores)
}