#' Convert a list into a GMT file
#'
#' This function converts a list into a GMT file
#' @param x list
#' @param gmtFile where the GMT file should be saved
#' @param needsNames logical, does the GMT file need fake names in
#' the first two positions of each row?
#' @title Convert a list into a GMT file
#' @author Giovanni d'Ario
#' @export
#' @examples
#' lst <- list(a = c(1,2), b = c(3,4,5))
#' gmt <- listToGmt(lst)
list_to_gmt <- function(x, gmtFile, needsNames=TRUE) {
    i <- 1
    if (needsNames) {
        nms <- paste0("S", i)
        x[[1]] <- c(nms, nms, x[[1]])
    }
    cat(file = gmtFile, paste(x[[1]], sep = "",
        collapse = "\t"), "\n")
    for(i in 2:length(x)) {
        if (needsNames) {
            nms <- paste0("S", i)
            x[[i]] <- c(nms, nms, x[[i]])
        }
        cat(file = gmtFile, paste(x[[i]], sep = "",
                collapse = "\t"),
            "\n", append = TRUE)
    }
}
