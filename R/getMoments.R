#' Calculate moments of a dataset
#' 
#' getMoments takes an input of a matrix of signatures, where the rows are features and the
#' columns are signatures.  It then calculates the mean, variance, and skewness of each row and applies 
#' a test for normality. 
#' 
#' @param mymat Numeric matrix of dimension features x samples
#' @param outpath Output path for the rds object with the calculated data. Default: ./moments.rds
#' @param rownames Optional character vector of character identifiers for the rows of the matrix. 
#' 
#' @returns A dataframe with row identifiers, the moments, and the normality test
#' @export
#' @importFrom stats var
#' @importFrom moments skewness
getMoments <- function(mymat, outpath="./moments.rds", rownames=c()){
  
  if (length(rownames) == 0){
    rownames <- row.names(mymat)
  } 
  if (is.null(rownames)){
    rownames <- seq_len(dim(mymat)[1])
  }
  
  mymeans <- rowMeans(mymat)
  
  myvar <- apply(mymat, 1, stats::var)
  
  myskew <- apply(mymat, 1, moments::skewness)
  myskew95 <- sapply(seq_len(dim(mymat)[1]), FUN=function(x) 
                moments::skewness(sort(mymat[x,])[round(dim(mymat)[2]*0.05): round(dim(mymat)[2]*0.95)]))
  
  #myshapiro <- apply(mymat, 1, shapiro.test)
  #normTestPval <- sapply(myshapiro, FUN=function(x) x$p.value)
  #normTestStat <- sapply(myshapiro, FUN=function(x) x$statistic)
  
  momentDF <- data.frame(rowID=rownames, 
                         mean=mymeans, 
                         var=myvar,
                         skew=myskew, 
                         skew95=myskew95)
                         #shapiroStat=normTestStat,
                         #shapiroPval=normTestPval)
  
  saveRDS(momentDF, file=outpath)
  
  return(momentDF)
}