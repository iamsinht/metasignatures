#library(dplyr)
#library(gridExtra)
#library(perturbKit)
#library(pracma)
#library(MASS)


# This function is used to sample a multivariate normal distribution
# with eigenvalues drawn from some distribution and compute var(cos).
# This can then be used to test the variance conjecture. 
# Add parameters to allow specification of the function from which
# eigenvalues are sampled.

#' Sample the variance of cosine similarity for vectors drawn from a multivarite normal
#' distribution with eigenvalues drawn from a Gamma distribution. 
#' 
#' @param ndim numeric, a vector of dimensions to sample. Default seq(100, 1000, 100)
#' @param iter Number of iterations to use for a particular dimension and eigenvalue distribution
#' @param distiter Number of eigenvalue distributions to run for each dimension
#' @param sPar Numeric, vector of gamma distribution hyperparameters for the eigenvalues. The eigenvalues
#' are sampled from a gamma distribution with shape parameter gamma(shape=sPar[1], rate=sPar[2]). and rate
#' parameter gamma(shape=sPar[3], rate=sPar[4])
#' 
#' @returns Dataframe with columns: dimension, theoretical variance based on eigenvalues, theoretical variance 
#' on observed eigenvalues, and empirical variance. 
#' @export
#' @importFrom stats rgamma
sampleVarCos <- function(ndim=seq(100, 1000, 100), iter=3, distiter=10, sPar=c(1, 2, 1, 2)){
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), theoryVar=numeric(mylen), theoryEmpVar=numeric(mylen), obsVar=numeric(mylen))
  k <- 1
  
  for (mydim in ndim){
    for (ii in seq_len(distiter)){
      print(sprintf("ii = %d", ii))
      a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
      b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
      
      for (jj in seq_len(iter)){
        sigmat <- diag(rgamma(mydim, shape=a, rate=b))
        lambda <- diag(sigmat)
        
        x <- t(MASS::mvrnorm(n=2000, mu=numeric(mydim), Sigma=sigmat))
        xpr <- stats::prcomp(t(x))
        
        xcos <- perturbKit::cosine(x, x)
        
        ovar <- stats::var(xcos[upper.tri(xcos)])
        tvar <- sum(lambda ** 2)/(sum(lambda) ** 2)
        
        tempvar <- sum(xpr$sdev^2 ** 2)/(sum(xpr$sdev^2))^2
        
        resdf[k,] <- c(mydim, tvar, tempvar, ovar)      
        k <- k + 1
      }
    }
  }
  
  return (resdf)
  
}
