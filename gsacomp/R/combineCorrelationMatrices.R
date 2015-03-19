createCorrelationMatrix <- function(size, rho) {
    if(rho > 1)
        stop("rho must be <= 1")
    M <- matrix(rho, nrow = size, ncol = size)
    diag(M) <- 1
    return(M)
}

##' Create a block correlation matrix
##'
##' This function creates a block correlation matrix where each
##' block can have a different size. Within each block, the
##' correlation is the same and equals rho. The correlation
##' among different blocks is zero.
##' @title Create a diagonal block correlation matrix
##' @param sizes numeric vector, the sizes of the individual
##' correlation matrices
##' @param rho numeric, the correlation coefficient
##' @return a diagonal block matrix
##' @author Giovanni d'Ario
##' @export
createBlockCorrelationMatrix <- function(sizes, rho) {
    require(Matrix)
    matrixList <- lapply(sizes, createCorrelationMatrix, rho = rho)
    result <- do.call(bdiag, matrixList)
    return(as.matrix(result))
}
