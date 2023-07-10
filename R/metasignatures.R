#' Extract metasignatures
#' getMetasigs takes an input of a matrix of signatures, where the rows are features and the
#' columns are signatures.  It generates metasignatures by averaging k signatures together, where
#' k is the given parameter (default 1). 
#' 
#' @param mymat Numeric matrix of dimension features x samples
#' @param k Numeric, number of signatures to sample and average without replacement
#' @param return2 logical: if true, it returns a list of two disjoint metasignatures of size k. 
#' @param returnk logical: if true and return2 is false, returns as many disjoint 
#' metasignatures of size k as possible. If return2 and returnk are both false, returns one metasignature
#' 
#' @returns some number of metasignatures depending on parameters specified
#' @export
#' @importFrom cmapR GCT
getMetasigs <- function(mymat, k=1, return2=1, returnk=0){

  if (k==1){
    if (return2 == 0){
      return(mymat[, sample(dim(mymat)[2], 1)])
    } else {
      ix <- sample(dim(mymat)[2], 2)
      return(cbind(mymat[, ix[1]], mymat[, ix[2]]))
    }
  }
  
  if (return2 == 0){
    ix <- sample(dim(mymat)[2], k)
    return(rowMeans(mymat[, ix]))
  } else if (returnk==1) {
    
    ix <- sample(dim(mymat)[2], floor(dim(mymat)[2]/k) * k)
    return(sapply(seq(floor(dim(mymat)[2]/k)), FUN=function(x) rowMeans(mymat[, ix[(k*(x-1)+1):(k*x)]])))

  } else {
    if (2*k <= dim(mymat)[2]){
      ix <- sample(dim(mymat)[2], 2*k)
      return(sapply(seq(2), FUN=function(x) rowMeans(mymat[, ix[(k*(x-1)+1):(k*x)]])))
      #return(list(rowMeans(mymat[, ix[1:n]]), rowMeans(mymat[, ix[(n+1):(2*n)]])))
    } else {
      return(c(0,0))
    }
  }
  
}


#' getMetaCorr
#' For an input matrix, compute the Pearson's correlation of metasignatures drawn from that matrix
#' as a function of the metasignature size k. 
#' @param mymat Numeric matrix of dimension features x samples
#' @param kmax Numeric, maximum value to use for metasignature size k. 
#' @param iter Numeric, the number of pairs of disjoint metasignatures to sample for each value of 
#' metasignature size k. Default = 100. 
#' @param useall Logical, whether to sample a range of values for kmax or use all values. Default FALSE, 
#' deprecated. 
#' 
#' @returns A list of (1) cormat - a matrix with dimension kmax x iter, where each entry is the Pearson's
#' correlation of a metasignature of size k, 
#' (2) cordf - a data frame with columns metasize, mean correlation, sd correlation summarizing the 
#' correlation of metasignatures generated from the input matrix. 
#' @export
getMetaCorr <- function(mymat, kmax, iter=100, useall=FALSE){
  kmax <- min(kmax, floor(dim(mymat)[2]/2))
  
  # Revisit:
  # if (kmax > 50){
  #   seqvals <- c(seq(10), seq(15, kmax, 5))
  #   seqvals <- seqvals[seqvals <= kmax]
  # } else {
  #   seqvals <- seq(kmax)
  # }
  
  # if (!useall){
  #   seqvals <- c(seq(10), seq(15, kmax, 5))
  #   seqvals <- seqvals[seqvals <= kmax]
  # } else {
  #   seqvals <- seq(1, kmax)
  # }
  seqvals <- seq(1, kmax, 1)
  
  cormat <- matrix(numeric(iter*length(seqvals)), nrow=length(seqvals), dimnames=list(seqvals))
  
  for (jj in seq(length(seqvals))){
    ii <- seqvals[jj]
    a <- lapply(seq(iter), FUN=function(x) getMetasigs(mymat, k=ii, return2=1))
    cormat[jj,] <- sapply(a, FUN=function(x) stats::cor(x)[2])
    #cormat[ii,] <- sapply(seq(iter), FUN=function(x) cor(a[[1,x]], a[[2,x]]))
  }
  
  cordf <- data.frame(metasize=seqvals, meancor=rowMeans(cormat), sdcor=apply(cormat, 1, stats::sd))
  return(list(cormat=cormat, cordf=cordf))
}


#' bootstrapMaxMetaCorr
#' For an input matrix, bootstraps the correlation of maximal disjoint metasignatures, i.e.
#' those of size floor(dim(mymat)[2]/2). 
#' @param mymat Numeric matrix of dimension features x samples
#' @param iter Numeric, the number of pairs of disjoint metasignatures to sample
#' 
#' @returns Numeric vector of length iter corresponding to the Pearson correlations of maximal disjoint
#' metasignatures sampled from the matrix. 
#' @export
bootstrapMaxMetaCorr <- function(mymat, iter=100){
  kmax <- floor(dim(mymat)[2]/2)
  
  mycors <- sapply(seq(iter), FUN=function(x) stats::cor(getMetasigs(mymat, k=kmax, return2=1))[2])
  return(mycors)
}



# This function is inordinately slow (because of wtcs)
#' getMetaSim
#' For an input matrix, compute the similarity of metasignatures drawn from that matrix
#' as a function of the metasignature size k. Analogous to getMetaCorr, but allows Pearson, 
#' Spearman, cosine similarity, and L1000's weighted connectivity score (wtcs).
#' @param mymat Numeric matrix of dimension features x samples
#' @param kmax Numeric, maximum value to use for metasignature size k. 
#' @param iter Numeric, number of pairs of disjoint metasignatures to sample for each value of 
#' metasignature size k. Default = 100. 
#' @param metric One of "pearson", "spearman", "wtcs", or "cosine". 
#' 
#' @returns List of: simmat - similarity matrix of dimension kmax x iter with the similarities 
#' of the metasignatures; cordf - data frame summarizing the mean and standard deviation of the
#' similarities; metric - string indicating which metric was used. 
#' @export
#' @importFrom coop cosine
getMetaSim <- function(mymat, kmax, iter=100, metric="pearson"){
  
  metric <- match.arg(metric, c("pearson", "spearman", "wtcs", "cosine"))
  kmax <- min(kmax, floor(dim(mymat)[2]/2))
  
  seqvals <- seq(kmax)
  
  # if (kmax > 50){
  #   seqvals <- seq(5, kmax, 5)
  # } else {
  #   seqvals <- seq(kmax)
  # }
  
  simmat <- matrix(numeric(length(seqvals)*iter), nrow=length(seqvals))
  
  myfunc <- switch(metric, 
                   pearson = function(x) stats::cor(x)[2], 
                   spearman = function(x) stats::cor(x, method="spearman")[2], 
                   wtcs = function(x) perturbKit::calcSimBlock(cmapR::GCT(mat=x, cid=as.character(seq(dim(x)[2]))), 
                                                               cmapR::GCT(mat=x, cid=as.character(seq(dim(x)[2]))), metric=metric)[2], 
                   cosine = function(x) coop::cosine(x)[2]) 
                   #wtcsmat = function(x) compute_cs_mats(x, x)[2])  # <- wtf is this
  
  for (ii in seq_along(seqvals)){
    a <- lapply(seq(iter), FUN=function(x) getMetasigs(mymat, k=seqvals[ii], return2=1))
    simmat[ii,] <- sapply(a, FUN=myfunc)
  }
  
  cordf <- data.frame(metasize=seqvals, meanSim=rowMeans(simmat), sdSim=apply(simmat, 1, stats::sd))
  rownames(simmat) <- seqvals
  return(list(simmat=simmat, cordf=cordf, metric=metric))
}


#' getMetaSimDs
#' Apply getMetaSim to each group within a dataset. 
#' @param mymat Numeric matrix of dimension features x samples
#' @param groupings Vector of class labels, either character or numeric. For each class, metacorrelation 
#' is computed.
#' @param kmax Numeric, maximum value to use for metasignature size k (default 1000). For any group with 
#' fewer than 2*kmax signatures, computes up to floor(groupsize/2). 
#' @param iter Numeric, number of pairs of disjoint metasignatures to sample for each value of metasize k. 
#' Default  = 100.
#' @param metric One of "pearson" (default), "spearman", "wtcs", or "cosine". 
#' 
#' @export
getMetaSimDs <- function(mymat, groupings, kmax=1000, iter=100, metric="pearson"){
  
  if (dim(mymat)[2] != length(groupings)){
    return("Error: length of groupings does not match number of columns of input matrix.")
  }
  
  ret <- list()
  mygrps <- names(table(groupings)[table(groupings) >= 6])
  
  for (ii in seq_along(mygrps)){
    if (ii %% 50 == 0){
      print(sprintf("%d/%d", ii, length(mygrps)))
    }
    agrp <- mygrps[ii]
    ix <- which(groupings == agrp)
    
    metaStr <- getMetaSim(mymat[, ix], kmax=kmax, iter=iter, metric=metric)
    ret <- c(ret, list(metaStr))
  }
  
  names(ret) <- mygrps
  
  return(ret)
}


#' ecdfPointwise
#' Takes a numeric vector and computes an ecdf on user-input values of x. This is particularly 
#' useful for overlaying multiple ECDFs at the same points. 
#' @param mydat numeric vector to compute an ECDF
#' @param breaks vector of points at which to compute the ECDF.
#' 
#' @returns A dataframe with columns x - the domain, the breaks; and y - the fraction of input 
#' values less than x. 
#' @export
ecdfPointwise <- function(mydat, breaks){
  h <- graphics::hist(mydat, breaks=breaks, plot=FALSE)
  y0 <- mean(mydat < breaks[1])
  return(data.frame(x=breaks, y=c(y0, cumsum(h$counts)/length(mydat))))
}
