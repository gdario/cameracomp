#' Read the output of a Camera simulation
#' 
#' Read and put together the results of a Camera simulation
#' 
#' This function reads the RData file produced by 
#' \code{run_simulation} and assembles a data frame. 
#' @param A string containing the full path to the output of 
#' \code{run_simulation}
#' @export
#' @author Giovanni d'Ario
read_camera_results <- function(resultsFile, outFile=NULL) {
  rf <- load(resultsFile)
  assign("out", get(rf))
  out <- lapply(out, function(x) {
    x$gene_set <- as.numeric(sub("set_", "", rownames(x)))
    rownames(x) <- NULL
    x
    })
  camera_results <- do.call(rbind, out)
  outDir <- dirname(resultsFile)
  outFile <- file.path(outDir, "camera_results.RData")
  save(camera_results, file = outFile)
  camera_results
}