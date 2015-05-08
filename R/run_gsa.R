##' Create the input for gsa.py
##'
##' This function takes the output of \code{createInputList} and turns it
##' into the input for GSA
##' @title Create the input for GSA
##' @param inputList list, the output of \code{createInputList}.
##' @return a data frame containing the output of limma on the simulated
##' dataset
##' @author Giovanni d'Ario
create_gsa_input <- function(inputList, adj) {
    require(limma)
    design <- model.matrix(~ inputList$y)
    fit <- lmFit(inputList$X, design)
    fit <- eBayes(fit)
    ttab <- topTable(fit, coef = 2, n = Inf)
    out <- data.frame(gene_id = rownames(ttab),
                      fc = ttab$logFC,
                      pval = ttab$adj.P.Val)
    if (!adj)
        out$pval <- ttab$P.Value
    return(out)
}

##' Run gsa.py on a simulated dataset
##'
##' This function creates a gsa.py input file and a gmt file
##' starting from the output of \code{createInputList}
##' @title Run gsa.py on a simulated dataset
##' @param inputList 
##' @param resultDir 
##' @param prefix 
##' @param adj 
##' @return no value returned, but as a side effect gsa.py is run
##' @author Giovanni d'Ario
run_gsa <- function(inputList, resultDir, prefix, adj=TRUE) {
    inputfile <- create_gsa_input(inputList, adj = adj)
    write.table(inputfile, file = file.path(resultDir, "infile.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    gmt <- list_to_gmt(x = inputList$index,
                     gmtFile = file.path(resultDir, "gmtfile.gmt"))
    gmtfile <- file.path(resultDir, "gmtfile.gmt")
    infile <- file.path(resultDir, "infile.txt")
    prefix <- file.path(resultDir, prefix)
    system(paste("gsa.py -i", infile, "-g", gmtfile, "--minsize 3",
                 "--common 3 --maxsize 600 -p", prefix))
}
