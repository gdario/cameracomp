##' Create the input list for Camera and GSA
##' 
##' @title Create a simulated input for Camera and GSA
##' @param sizes integer vector, the sizes of the 
##' @param rho double, the within-gene-set correlation coefficient
##' @param nsamples integer, the total number of samples.
##' @return a list containing the expression matrix X, the class
##' vector y, and the list of gene sets `genesets`
##' @author Giovanni d'Ario
create_input_list <- function(sizes, rho, nsamples) {
    X <- create_gene_expression_matrix(sizes = sizes, rho = rho,
                                    nsamples = nsamples)
    ## rownames(X) <- paste0("g", 1:nrow(X))
    y <- gl(n = 2, k = nsamples/2, labels = c("A", "B"))

    genesets <- data.frame(gene = 1:nrow(X),
                           set = rep(1:length(sizes), sizes))
    genesets <- split(genesets$gene, genesets$set)
    ## genesets <- lapply(genesets, as.character)

    return(list(X = X, y = y, index = genesets))
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


