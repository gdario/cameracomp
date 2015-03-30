#' Create a correlation matrix with a given size and correlation
#' 
#' This function takes in input an integer \code{n} indicating the
#' number of rows (and columns, since the matrix is symmetrical) 
#' and a real number indicating the pairwise correlation. It then
#' returns a \code{n x n} matrix with ones on the diagonal and 
#' \code{rho} off the diagonal
#' @title Create a correlation matrix
#' @param size. Integer, the number of rows (columns) of the matrix.
#' @param rho. Real, the pairwise correlation coefficient.
#' @author Giovanni d'Ario
create_correlation_matrix <- function(size, rho) { 
    if (rho < 0 | rho > 1)
        stop("rho must be a number beetween 0 and 1")
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
create_block_matrix <- function(sizes, rho) {
    require(Matrix)
    matrixList <- lapply(sizes, create_correlation_matrix, rho = rho)
    result <- do.call(bdiag, matrixList)
    return(as.matrix(result))
}
