#' Create a set of differentially expressed (and not) gene sets
#' 
#' Create a gene expression matrix containing a number of 
#' differentially expressed and of non-differentially expressed
#' gene sets, plus a number of non-differentially expressed genes
#' not included in any gene set.
#' 
#' This function produces a gene expression matrix of \code{n_genes}
#' genes part of which are members of a gene set. More precisely
#' \code{n_de_sets} gene sets contain a proportion \code{p_de} of
#' differentially expressed gene (where by differentially expressed
#' we mean that their fold change is \code{logfc}). In addition, the
#' expression matrix contains \code{n_nde_sets} gene sets each 
#' containing exclusively non-differentially expressed genes. These
#' sets exist to measure the impact of the gene set size on the 
#' number of false positives. The size of each gene set can be 
#' chosen by modifying the parameters \code{size_de_sets} and 
#' \code{size_nde_sets} respectively. 
#' 
#' The expression matrix
#' is divided into two balanced groups, for a total of 
#' \code{n_samples} samples. Each gene set (regardless of being
#' differentially expressed) has an intra-set correlation of 
#' \code{rho}. The genes that do not belong to any gene set have 
#' zero correlation.
#' 
#' @param n_de_sets Integer. Number of differentially expressed gene
#' sets.
#' @param size_de_sets Integer vector. The sizes of the 
#' \code{n_de_sets} differentially expressed gene sets.
#' @param n_genes Integer. The total number of genes in the
#' expression matrix. Default is 10000.
#' @param n_nde_sets Integer, the number of non-differentially 
#' expressed gene sets.
#' @param size_nde_sets Integer vector. The sizes of the 
#' \code{n_nde_sets} non-differentially expressed genes.
#' @param n_samples Integer, the total number of 
#' subjects. The final expression matrix will have size
#'  \code{n_genes} * \code{n_samples}.
#' @param p_de Real. The proportion of genes in the gene set that
#' are differentially expressed. Default is 50\%.
#' @param rho Real. The intra-set correlation.
#' @param logfc Real. The average fold change in the differentially
#' expressed genes.
#' @export
#' @author <Giovanni d'Ario giovanni.dario@@novartis.com>
#' @examples 
#' input_list <- create_de_gene_sets(n_de_sets = 5, rho = 0.1)
create_de_gene_sets <- function(n_de_sets=5L,
                          size_de_sets=seq(from=20, 
                                           to = 500,
                                           length.out = n_de_sets),
                          n_nde_sets=n_de_sets,
                          size_nde_sets=size_de_sets,
                          n_genes=10000L,
                          n_samples=40L,
                          p_de=0.5,
                          rho=.5, 
                          logfc=1.0) {
  
  ## The same covariance matrix is used for the differentially
  ## expressed and for the non differentially expressed gene 
  ## sets
  S <- create_block_matrix(sizes = size_de_sets, rho = rho)
  
  ## Matrix of differential expressed genes
  group1 <- MASS::mvrnorm(n = .5 * n_samples,
                          mu = rep(0, sum(size_de_sets)),
                          Sigma = S)
  group1 <- t(group1)
  
  group2 <- lapply(seq_len(n_de_sets), function(i) {
    s <- create_correlation_matrix(size = size_de_sets[i], rho = rho)
    ## p_de * setsize genes are differentially expressed
    nde <- round(p_de * size_de_sets[i])
    m <- c(rep(logfc, nde), rep(0, size_de_sets[i] - nde))
    g <- MASS::mvrnorm(n = .5 * n_samples, mu = m, Sigma = s)
    t(g)
  })
  group2 <- do.call(rbind, group2)
  
  ## Expression matrix for the differentially expressed gene sets
  Mde <- cbind(group1, group2)
  
  ## Create the non-differentially expressed gene sets
  Snde <- create_block_matrix(sizes = size_nde_sets, rho = rho)
  Mnde <- MASS::mvrnorm(n = n_samples,
                        mu = rep(0, sum(size_nde_sets)),
                        Sigma = Snde)
  Mnde <- t(Mnde)
  
  ## Matrix of non-differentially expressed genes
  n_others <- n_genes - (nrow(Mde) + nrow(Mnde))
  
  other_genes <- matrix(rnorm(n_samples * n_others),
                        nrow = n_others,
                        ncol = n_samples)
  M <- rbind(Mde, Mnde, other_genes)
  gene_set <- rep(paste0("set_", seq_len(n_de_sets + n_nde_sets)),
                  c(size_de_sets, size_nde_sets))
  gene_set_df <- data.frame(gene = seq_along(gene_set),
                            gene_set = gene_set)
  index <- split(gene_set_df$gene, gene_set_df$gene_set)
  
  y = gl(n = 2, k = .5 * n_samples, labels = c("A", "B"))
  
  out <- list(X = M, y = y, index=index)
  
  out
}