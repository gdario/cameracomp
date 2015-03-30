##' Run a simulation on Camera and gsa.py
##'
##' This function runs a user-specified number of simulations creating
##' artificial datasets with no differentially expressed genes, but
##' a block-diagonal correlation matrix. The output of camera is
##' saved into an RData object for further inspection. The output
##' of gsa.py is saved in the \code{resultDir} in the form of
##' multiple files, one for each run of the simulation.
##' @title Run a simulation using camera and gsa.py on synthetic data
##' @param sizes integer vector containin the sizes of the gene ses
##' @param rho float, the within-gene-set correlation coefficient
##' @param nsamples integer, the total number of samples in the
##' data set
##' @param resultDir string, the path to the directory where the
##' results will be saved.
##' @param B integer, the number of runs in the simulation
##' @return a list containing the results of the camera simulation
##' (which are also saved as an RData object) and a set of files
##' from gsa.py
##' @author Giovanni d'Ario
##' @export
run_simulation <- function(sizes, rho, nsamples, resultDir=NULL, B) {

    if (is.null(resultDir))
        resultDir <- paste0("results_", rho, "_b_", B)
    if (!file.exists(resultDir))
        dir.create(resultDir)
    
    out <- lapply(1:B, function(i) {
        message(paste("B =", i))
        inputList <- create_input_list(sizes = sizes,
                                     rho = rho,
                                     nsamples = nsamples)
        res <- run_camera(inputList)
        run_gsa(inputList, adj = TRUE, resultDir = resultDir, prefix = i)
        return(res)
    })
    message("Simulation complete")
    save(out, file = file.path(resultDir, "/results_simulation.RData"))
    return(invisible(out))
}

clean_all <- function(resultDir) {
    txtFiles <- file.path(resultDir, "*.txt")
    gmtFiles <- file.path(resultDir, "*.gmt")
    rdataFiles <- file.path(resultDir, "*.RData")
    system(paste("rm -f", txtFiles))
    system(paste("rm -f", gmtFiles))
    system(paste("rm -f", rdataFiles))
}
