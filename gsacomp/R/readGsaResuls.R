##' Read the output files of GSA and combine them in one object
##'
##' This function reads the n.txt files in the specified directory,
##' where n is an integer ranging from 1 to B, where B is the number
##' of simulations specified in \code{runSimulations}. The files are
##' then organized in one single data frame
##' @title Read the output files of GSA and combine them in one object
##' @param dirName string, the full path to the directory containing
##' the output of the GSA simulation
##' @param outFile string, the name of the output file. This file
##' will be automatically saved in the \code{dirName} folder.
##' @param nFiles integer, the number of files. If NULL it is
##' extracted from the output folder name, which is in the form
##' "results_rho_b_nFiles".
##' @return a data frame that combines all the results of the GSA
##' simulation. As a side effect the output will be saved in the
##' \code{dirName} folder under the name specified in \code{outFile}
##' as an RData file.
##' @author Giovanni d'Ario
readGsaResults <- function(dirName,
                           outFile=NULL,
                           nFiles=NULL) {

    if (is.null(outFile))
        stop("Please speficy a name and a location for the output file")
    
    if (is.null(nFiles))
        nFiles <- as.numeric(sub("^.*b_([0-9]+)/?$", "\\1", dirName))

    ## Extract the names of the various columns (the second row
    ## contains the column type and we discard it
    header <- names(read.delim(file.path(dirName, "1.txt"), nrows = 1))

    ## Read a file skipping the irrelevant information
    readFiles <- function(x, h) {
        tmp <- read.delim(x, skip = 1, stringsAsFactors = FALSE)
        names(tmp) <- h
        return(tmp)
    }
    
    gsaResults <- lapply(1:nFiles, function(i) {
        fileName <- file.path(dirName, paste(i, ".txt", sep = ""))
        tmp <- readFiles(x = fileName, h = header)
    })
    resultsGsa <- do.call('rbind', gsaResults)
    save(resultsGsa, file = file.path(dirName, outFile))
    return(resultsGsa)
}

##' Select one particular test from the GSA results
##'
##' This function extracts one particular test from the GSA
##' results and the associated log10 pavalues
##' @title Select one particular test from the GSA results
##' @param resultsFile path to the results file for GSA
##' @param column the name of the column(s) to be extracted.
##' Dafaults to 'fc_tTest'
##' @return a data frame containing the selected columns and, as
##' a side effect, an RData file containing the same data
##' @author Giovanni d'Ario
extractResultsFromGsa <- function(resultsFile, column="fc_tTest") {
    rf <- load(resultsFile)
    ## The name of the results file can change. Assign it
    ## to a variable called 'res'
    assign("res", get(rf))
    selectedColumns <- grep(column, names(res), value = TRUE)
    selectedColumns <- c("name", "gene_set_size", selectedColumns)
    selectedGsaResults <- res[selectedColumns]

    ## Put the p-values on the linear scale
    pvalNegLog <- selectedGsaResults[[grep("pval_neglog10",
                                           names(selectedGsaResults))]]
    adjPvalNegLog <- selectedGsaResults[[grep("pval_adj_neglog10",
                            names(selectedGsaResults))]]

    selectedGsaResults <- within(selectedGsaResults, {
        NGenes <- gene_set_size
        PValue <- 10 ^ (-pvalNegLog)
        FDR <- 10 ^ (-adjPvalNegLog)
    })

    selectedGsaResults <- subset(selectedGsaResults,
                                 select = c("NGenes", "PValue", "FDR"))
    
    outFile <- sub("\\.RData", "_selection.RData", resultsFile)
    save(selectedGsaResults, file = outFile)
    return(selectedGsaResults)
}

##' Create a data frame to compare Camera and GSA results
##'
##' This function combines the output from GSA and Camera in order to
##' allow the visualization of the results.
##' @title Create a data frame to compare Camera and GSA results
##' @param cameraFile string, path to the file containing the Camera
##' output
##' @param gsaFile string, path to the file containin the GSA output
##' @return A data frame containing the gene set size, the p-values,
##' the adjusted p-value and the method
##' @author Giovanni d'Ario
compareResults <- function(cameraFile, gsaFile) {
    cameraResults <- load(cameraFile)
    assign("cameraResults", get(cameraResults))
    gsaResults <- load(gsaFile)
    assign("gsaResults", get(gsaResults))

    cameraResults <- do.call('rbind', cameraResults)
    cameraResults <- subset(cameraResults,
                            select = c("NGenes", "PValue", "FDR"))

    cameraResults$method <- "camera"
    gsaResults$method <- "gsa_py"

    out <- rbind(cameraResults, gsaResults)

    return(out)
}
