#' Create a differentially expressed gene set and its complement
#' 
#' Create a gene expression matrix of non regulated
#' genes a subset of which is (in part) differentially expressed.
#' 
#' This function produces a gene expression matrix of \code{ngenes}
#' genes of in which a gene set of size \code{setsize} is 
#' differentially expressed. The gene set has an intra-set 
#' correlation of \code{rho} and a fraction \code{p_de} of 
#' its genes is differentially expressed, with an average
#' fold change of \code{logfc}. The genes in the complement
#' are not differentially expressed and have zero correlation.
#' 
#' @param ngenes Integer. The total number of genes in the
#' expression matrix. Default is 10000.
#' @param nsamples Integer, the total number of 
#' subjects. The final expression matrix will have size
#'  \code{ngenes} * \code{nsamples}.
#' @param sizes Integer. The size of the gene set. Default is 50.
#' @param p_de Real. The proportion of genes in the gene set that
#' are differentially expressed. Default is 50%.
#' @param rho Real. The intra-set correlation.
#' @param logfc Real. The average fold change in the differentially
#' expressed genes.
create_de_set <- function(ngenes=10000L,
                          nsamples=40L,
                          setsize=50L,
                          sizes=c(setsize, ngenes-setsize),
                          p_de=0.5,
                          rho=.5, 
                          logfc=1.0) {
  
  S <- create_correlation_matrix(size = setsize, rho = rho)
  
  ## Matrix of differential expressed genes
  group1 <- MASS::mvrnorm(n = .5 * nsamples,
                          mu = rep(0, setsize),
                          Sigma = S)
  ## p_de * setsize genes are differentially expressed
  nde <- round(p_de * setsize)
  m <- c(rep(logfc, nde), rep(0, setsize - nde))
  
  group2 <- MASS::mvrnorm(n = .5 * nsamples,
                          mu = m,
                          Sigma = S)
  Mde <- cbind(t(group1), t(group2))
  
  ## Matrix of non-differentially expressed genes
  Mnde <- matrix(rnorm(nsamples * (ngenes - setsize)),
                nrow = ngenes - setsize, 
                ncol = nsamples)
  M <- rbind(Mde, Mnde)
  out <- list(M = M, 
              setsize = setsize,
              rho = rho, 
              gene_set <- factor(c(rep("gene_set", setsize),
                                   rep("complement", ngenes - setsize))),
              cl = gl(n = 2, 
                      k = .5 * nsamples, 
                      labels = c("A", "B")))
  out
}