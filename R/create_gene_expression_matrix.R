##' Generate a simulated expression matrix
##' 
##' Create a simulated expression matrix with block correlations
##'
##' This function creates a fake expression matrix where the genes
##' are correlated according to the correlation matrix specified
##' by the block sizes and the correlation coefficient \code{rho}.
##' @title Simulate an expressin matrix a with given covariance structure
##' @param sizes integer vector, the sizes of each gene set
##' @param rho real, the within-gene-set correlation
##' @param nsamples integer, the number of samples
##' @return a simulated expression matrix
##' @author Giovanni d'Ario
##' @export
create_gene_expression_matrix <- function(sizes, rho, nsamples) {
    S <- create_block_matrix(sizes = sizes, rho = rho)
    X <- MASS::mvrnorm(n = nsamples, mu = rep(0, nrow(S)), Sigma = S)
    return(t(X))
}
