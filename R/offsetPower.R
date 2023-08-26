#' offset Power 
#' 
#' This function is used to simulate adding an offset to L1000 data and assessing
#' the power to detect same-compound signatures. It is improvised for time and
#' should not be exported with the package. 
#' 
#' @param cellline The name of the cell line, e.g. "A375", "MCF7"
#' @param datapath The path to the directory containing the L1000 gctx files
#' @param metapath The path to the directory containing L1000 metadata
#' @param outpath The full RDS path to save the objects

offsetPower <- function(cellline="A375", datapath, metapath, outpath="./offsetpower.rds", iter=100, l1kmeta=c()){
  
  # Replace this as soon as humanly possible with a reference to metric learning.
  # Or probably better yet is migrate the functionality to perturbKit. 
  #source("~/Work/code/metriclearning/learningFunctions.R") # for innerProductGroups
  #source("~/Work/code/metriclearning/paperscripts/paperFigFuncs.R") # for balancedSample
  
  if (length(l1kmeta) == 0){
    l1kmeta <- CMAPToolkit::read_l1k_meta(datapath, version=2020)
  } else {
    print("L1k Metadata supplied, skipping reading")
  }
  
  resdf <- data.frame()
  deltadf <- data.frame()
  
  
  cellsiginfo <- l1kmeta$siginfo[l1kmeta$siginfo$cell_id == cellline & l1kmeta$siginfo$pert_type == "trt_cp", ]
  ds <-  cmapR::parse_gctx(perturbKit::get_level5_ds(datapath), cid=cellsiginfo$sig_id, rid=l1kmeta$landmarks$pr_gene_id)
  
  for (ii in seq(iter)){
    print(ii)
    
    # Pick a unit vector
    myvec <- rnorm(978)
    myvec <- myvec/sqrt(sum(myvec^2))
    
    # Determine the variance of the dataset along the unit vector
    x <- myvec %*% ds@mat
    xsd <- sd(x)
    
    deltadf <- rbind(deltadf, data.frame(cellid = cellline, dsMean = mean(x)/xsd))
    
    # Center
    matCtr <- ds@mat - mean(x) * matrix(rep(myvec, dim(ds@mat)[2]), ncol=dim(ds@mat)[2])
    
    for (delta in c(0, exp(seq(log(0.1), log(20), log(200)/12)))){
      matShift <- matCtr + delta * xsd * matrix(rep(myvec, dim(ds@mat)[2]), ncol=dim(ds@mat)[2])
      
      iGrps <- innerProductGroups("pearson", t(matShift), cellsiginfo$pert_iname, compact = 1)
      yr <- rankVectors(unlist(balancedSample(iGrps$same, k=100)), iGrps$diff)
      
      yh <- (unlist(balancedSample(iGrps$same)) - mean(iGrps$diff))/sd(iGrps$diff)
      
      resdf <- rbind(resdf, data.frame(cellid = cellline, ii = ii, delta = delta, 
                                       pLT01 = mean(yr < 0.01), 
                                       fdrLT10 = mean(p.adjust(yr, "fdr") < 0.1), 
                                       fdrLT25 = mean(p.adjust(yr, "fdr") < 0.25), 
                                       SNR2 = mean(yh > 2)))
    }
  }
  
  ret <- list(resdf=resdf, deltadf=deltadf, cellline=cellline)
  
  saveRDS(ret, file=outpath)
  
  return(ret)
}




## This is a huge hack, need to add to perturbKit:
## InnerProductGroups - grabbed explicitly from metriclearning/learningFunctions.R
innerProductGroups <- function(model, mat1, classes, compact=0){
  if (compact){
    # Use compact if the matrix is too large to effectively compute
    # Compact samples from the space of unlike similarities rather than computing the entire matrix
    
    # This could likely be accomplished more elegantly with dplyr
    samesim <- sapply(names(table(classes)[table(classes)>1]), FUN=function(x){
      ix <- which(classes == x)
      a <- innerProduct(model, mat1[ix, ], mat1[ix, ])
      a[upper.tri(a)]
    })
    
    jx <- sample(dim(mat1)[1], min(1000, dim(mat1)[1]))
    diffsim <- innerProduct(model, mat1[jx, ], mat1[jx, ])
    for (mygrp in names(table(classes)[table(classes) > 1])){
      ix <- which(jx %in% which(classes == mygrp))
      if (length(ix) > 1){
        diffsim[ix, ix] <- NA
      }
    }
    diffsim <- diffsim[upper.tri(diffsim)]
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
    
  } else {
    sims <- innerProduct(model, mat1, mat1)
    
    samesim <- sapply(names(table(classes)[table(classes) > 1]), FUN=function(x) {
      a <- sims[which(classes == x), which(classes == x)]
      a[upper.tri(a)]
    })
    
    diffsim <- sims
    for (mygroup in classes){
      ix <- which(classes == mygroup)
      diffsim[ix,ix] <- NA
    }
    
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
  }
}

## InnerProduct - grabbed explicitly from metriclearning/learningFunctions.R
innerProduct <- function(model, mat1, mat2){
  
  if (is.character(model)){
    if (model == "cosine"){
      return(cosine(mat1, mat2))
    }
    if (tolower(model) == "pearson"){
      return(cor(t(mat1), t(mat2), method="pearson"))
    }
    if (tolower(model) == "spearman"){
      return(cor(t(mat1), t(mat2), method="spearman"))
    }
  }
  # Check for integrity of variables?
  m1 <- model(torch_tensor(mat1, dtype=torch_float()))
  m2 <- model(torch_tensor(mat2, dtype=torch_float()))
  
  return(cosine(as.matrix(m1), as.matrix(m2)))
}

## balancedSample - taken from metriclearning/paperscripts/paperFigFuncs.R
balancedSample <- function(mylist, k=100){
  
  return(sapply(mylist, FUN=function(x) sample(x, min(length(x), k))))
  
}