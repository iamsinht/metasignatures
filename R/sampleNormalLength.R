#' Calculate the empirical Euclidean length of multivariate normal random vectors with eigenvalues sampled
#' from a gamma distribution.  
#'
#' @inheritParams sampleVarCos
#'
#' @returns Dataframe with columns: ndim - dimension of vector, meanLen - average length of vectors, sdLen - 
#' standard deviation of length, theoryLength - theoretical length predicted according to Central Limit Theorem, Smith et al.
#' @export
#' @importFrom stats sd
getNormalLength <- function(ndim=seq(100, 5000, 100), iter=10, distiter=10, sPar=c(1,2,1,2)){
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), meanLen=numeric(mylen), sdLen=numeric(mylen), theoryLength=numeric(mylen))
  k <- 1
  
  for (mydim in ndim){
    print(sprintf("mydim = %d", mydim))
    for (ii in seq_len(distiter)){
      print(sprintf("ii = %d", ii))
      a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
      b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
      
      for (jj in seq_len(iter)){
        sigmat <- diag(stats::rgamma(mydim, shape=a, rate=b))
        lambda <- diag(sigmat)
        
        x <- t(MASS::mvrnorm(n=100, mu=numeric(mydim), Sigma=sigmat))
        
        meanLen <- mean(sqrt(colSums(x^2)))
        sdLen <- stats::sd(sqrt(colSums(x^2)))
        theoryLen <- sqrt(sum(lambda))
        
        resdf[k,] <- c(mydim, meanLen, sdLen, theoryLen)
        k <- k + 1
      }
    }
  }
  
  return (resdf)
  
}

#' Calculate the empirical Euclidean length of multivariate normal random vectors with eigenvalues sampled
#' from a gamma distribution and non-zero mean sampled from a normal distribution
#'
#' @inheritParams sampleVarCos
#' @param meanSd Numeric, standard deviation of the normal distribution from which to sample means
#'
#' @returns Dataframe with columns: ndim - dimension of vector, meanLen - average length of vectors, sdLen - 
#' standard deviation of length, theoryLength - theoretical length predicted according to Central Limit Theorem, Smith et al.
#' @export
#' @importFrom stats sd
getGenNormalLength <- function(ndim=seq(100, 5000, 100), iter=10, distiter=10, sPar=c(1,2,1,2),
                               meanSd = 2){
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), meanLen=numeric(mylen), sdLen=numeric(mylen), 
                      theoryLength=numeric(mylen), tLenNoMean=numeric(mylen))
  k <- 1
  
  for (mydim in ndim){
    print(sprintf("mydim = %d", mydim))
    for (ii in seq_len(distiter)){
      print(sprintf("ii = %d", ii))
      a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
      b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
      
      for (jj in seq_len(iter)){
        sigmat <- diag(stats::rgamma(mydim, shape=a, rate=b))
        lambda <- diag(sigmat)
        
        myMeans <- stats::rnorm(n = mydim, mean = 0, sd = meanSd)
        x <- t(MASS::mvrnorm(n=100, mu=myMeans, Sigma=sigmat))
        
        meanLen <- mean(sqrt(colSums(x^2)))
        sdLen <- stats::sd(sqrt(colSums(x^2)))
        
        theoryLen <- sqrt(sum(myMeans^2 + lambda))
        tLenNoMean <- sqrt(sum(lambda))
        
        resdf[k,] <- c(mydim, meanLen, sdLen, theoryLen, tLenNoMean)
        k <- k + 1
      }
    }
  }
  
  return (resdf)
  
}
