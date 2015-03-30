##' Read a GMT file into a list
##'
##' This function reads a GMT file into an R list. At the moment
##' no type checking is performed, and all the elements in the
##' GMT file are assumed to be of type character. This may change
##' in the future, but it is important to be aware of this,
##' especially when working with ENTREZ ids.
##' @title Read a GMT file into a list
##' @param gmtFile full path to the GMT file
##' @return an R list where each component is a gene set
##' @author Giovanni d'Ario
gmt_to_list <- function(gmtFile) {
    gmt <- scan(gmtFile, sep = "\n", what = "character")
    gmt <- lapply(gmt, strsplit, split = "\t")
    unlist(gmt, recursive = FALSE)
}
