#library(dplyr)
#library(gridExtra)
#library(perturbKit)
#library(pracma)
# library(MASS)


# This function is used to sample a multivariate normal distribution
# with eigenvalues drawn from some distribution and compute var(cos).
# This can then be used to test the variance conjecture. 
# Add parameters to allow specification of the function from which
# eigenvalues are sampled.

#' Sample the variance of cosine similarity for vectors drawn from a multivariate normal
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


#' Sample the variance of cosine similarity for vectors drawn from a generalized multivariate normal
#' distribution with non-zero mean and with eigenvalues drawn from a Gamma distribution. 
#' 
#' @param ndim numeric, a vector of dimensions to sample. Default seq(100, 1000, 100)
#' @param iter Number of iterations to use for a particular dimension and eigenvalue distribution
#' @param distiter Number of eigenvalue distributions to run for each dimension
#' @param sPar Numeric, vector of gamma distribution hyperparameters for the eigenvalues. The eigenvalues
#' are sampled from a gamma distribution with shape parameter gamma(shape=sPar[1], rate=sPar[2]). and rate
#' parameter gamma(shape=sPar[3], rate=sPar[4])
#' @param meanSd Numeric, standard deviation of the normal distribution from which to sample means
#' @param nSamples Numeric, number of vectors to sample
#' 
#' @returns Dataframe with columns: dimension, theoretical variance based on eigenvalues, theoretical variance 
#' on observed eigenvalues, and empirical variance. 
#' @export
#' @importFrom stats rgamma
#' @importFrom stats rnorm
sampleGenVarCos <- function(ndim=seq(100, 1000, 100), iter=3, distiter=10, sPar=c(1, 2, 1, 2), 
                            meanSd=2, nSamples=2000){
  
  mylen <- length(ndim)*iter*distiter
  resdf <- data.frame(ndim=numeric(mylen), obsMean=numeric(mylen), theoryMean=numeric(mylen),
                      obsVar=numeric(mylen), theoryVar=numeric(mylen), tVarNoMean=numeric(mylen))
  k <- 1
  
  for (mydim in ndim){
    for (ii in seq_len(distiter)){
      print(sprintf("ii = %d", ii))
      a <- stats::rgamma(1, shape=sPar[1], rate=sPar[2])
      b <- stats::rgamma(1, shape=sPar[3], rate=sPar[4])
      
      for (jj in seq_len(iter)){
        sigmat <- diag(rgamma(mydim, shape=a, rate=b))
        lambda <- diag(sigmat)
        
        myMeans <- stats::rnorm(n = mydim, mean = 0, sd = meanSd)
        x <- t(MASS::mvrnorm(n=nSamples, mu=myMeans, Sigma=sigmat))
        
        #xpr <- stats::prcomp(t(x))
        
        xcos <- perturbKit::cosine(x, x)
        
        ovar <- stats::var(xcos[upper.tri(xcos)])
        tvar <- sum(lambda*(lambda + 2*(myMeans **2)))/(sum(myMeans**2 + lambda) ** 2)
        tvarNoMean <- sum(lambda ** 2)/(sum(lambda) ** 2)
        
        omean <- mean(xcos[upper.tri(xcos)])
        tmean <- sum(myMeans ** 2)/sum(myMeans**2 + lambda)
        #tempvar <- sum(xpr$sdev^2 ** 2)/(sum(xpr$sdev^2))^2
        
        resdf[k,] <- c(mydim, omean, tmean, ovar, tvar, tvarNoMean) 
        k <- k + 1
      }
    }
  }
  
  return (resdf)
  
}


#' Calculate the variance of cosine specifically for the two-dimensional case
#' 
#' @param myn - Number of vectors to sample at each ratio
#' @param iter - Number of different values of u to sample on [0, 0.5]
#' @param meanvec - Mean values to use (default (0,0)).
#' 
#' @returns dataframe with mean, variance, and approximation variance
#' @export
sampleVarCos2D <- function(myn=1000, iter=100, meanvec=c(0,0)){
  uvals <- seq(0.5/iter, 0.5, 0.5/iter)
  
  resdf <- data.frame(u=uvals, meanCos=numeric(length(uvals)), varCos=numeric(length(uvals)), theoryVar=numeric(length(uvals)))
  
  # S1, S2 are the eigenvalues of the covariance matrix, i.e. sigma_i^{2}
  for (myu in uvals){
    print(sprintf("iter = %d/%d", match(myu, uvals), iter))
    s1 <- myu
    s2 <- 1 - myu
    
    x <- t(MASS::mvrnorm(n=myn, mu=meanvec, Sigma=diag(c(s1, s2))))
    
    xcos <- perturbKit::cosine(x,x)
    
    meancos <- mean(xcos[upper.tri(xcos)])
    varcos <- stats::var(xcos[upper.tri(xcos)])
    theoryvar <- (s1^2 + s2^2)/((s1+s2)^2)

    resdf[match(myu, uvals), 2:4] <- c(meancos, varcos, theoryvar)
  }
  
  params <- c(myn=myn, iter=iter, meanvec=meanvec)
  
  return(list(resdf=resdf, params=params))
}


