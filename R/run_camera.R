##' Create the input list for Camera and GSA
##' 
##' @title Create a simulated input for Camera and GSA
##' @param sizes integer vector, the sizes of the 
##' @param rho double, the within-gene-set correlation coefficient
##' @param n_samples integer, the total number of samples.
##' @return a list containing the expression matrix X, the class
##' vector y, and the list of gene sets `genesets`
##' @author Giovanni d'Ario
create_input_list <- function(sizes,
                              de=FALSE,
                              n_de_sets=5L,
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
  
  if (de) {
    out <- create_de_gene_sets(n_de_sets = n_de_sets, 
                             n_nde_sets = n_nde_sets, 
                             n_genes = n_genes, 
                             p_de = p_de, 
                             size_de_sets = size_de_sets,
                             size_nde_sets = size_nde_sets, 
                             n_samples = n_samples, 
                             rho = rho, 
                             logfc = logfc)
  } else {
    X <- create_gene_expression_matrix(sizes = sizes, 
                                       rho = rho,
                                       n_samples = n_samples)
    
  ## rownames(X) <- paste0("g", 1:nrow(X))
  y <- gl(n = 2, k = n_samples/2, labels = c("A", "B"))
  
  genesets <- data.frame(gene = 1:nrow(X),
                         set = rep(1:length(sizes), sizes))
  genesets <- split(genesets$gene, genesets$set)
  ## genesets <- lapply(genesets, as.character)
  out <- list(X = X, y = y, index = genesets)
  }
  
  return(out)
}

##' Run camera on a simulated dataset
##'
##' This function runs camera on a simulated dataset produced by
##' \code{createCameraOutput}.
##' @title Run camera on a simulated dataset
##' @param inputList list, output of \code{createCameraInput}
##' @return a data frame containing the output of camera
##' @author Giovanni d'Ario
run_camera <- function(inputList) {
  design <- model.matrix(~ inputList$y)
  cameraOut <- limma::camera(inputList$X, inputList$index, design)
  return(cameraOut)
}