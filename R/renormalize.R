
#' renormData
#' This function takes a matrix and a subset of the columns, and renormalizes them using
#' one of several methods. The default is to center the matrix, so the mean is 0. It then
#' applies the transformation to a second matrix and returns the renormalized second matrix. 
#' In practice, the second matrix could be the same as the first matrix, in which case 
#' it returns the renormalized matrix. 
#' @param mat1 The matrix to center, of dimension features x samples
#' @param colix Numeric, the indices of columns to use (default is all). 
#' @param mat2 The matrix to return, of dimension features x samples
#' @param method Method to use (default = center). 
#' 
#' @returns A matrix, a renormalized mat2 per the criterion calculated on mat1. 
#' @export
renormData <- function(mat1, mat2, colix=c(), method="center"){
  
  if (dim(mat1)[1] != dim(mat2)[1]){
    stop("Error: dimensions of matrices do not match.")
  }

  if (length(colix) != 0){
    focusMat <- mat1[, colix]
  } else {
    focusMat <- mat1
  }
    
  if (method == "center"){
    shiftMat <- rowMeans(focusMat)
    retmat <- mat2 - matrix(rep(shiftMat, dim(mat2)[2]), ncol=dim(mat2)[2])
  }
  
  return(retmat)
}