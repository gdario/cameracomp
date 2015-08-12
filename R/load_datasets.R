#' Load the GSA and the Camera outputs at once
#' 
#' Load the Camera, the GSA t-test, KSw Dm and KSw Dp datasets from 
#' a given folder
#' 
#' This function receives the full path to a folder containing 
#' the output of a simulation, and loads the RData file with 
#' the results from Camera and from GSA's t-test, KSw Dm and Dp.
#' @return a list with one component for each of the tests 
#' mentioned above.
#' 
#' @param datadir A string containing the full path to a folder
#' storing the output of a simulation
#' @author Giovanni d'Ario <giovanni.dario@@novartis.com>
#' @export
load_datasets <- function(datadir) {
  load(file.path(datadir, "camera_results.RData"))
  load(file.path(datadir, "gsa_results.RData"))
  load(file.path(datadir, "gsa_results_KSw_Dm.RData"))
  load(file.path(datadir, "gsa_results_KSw_Dp.RData"))
  datasets <- list(camera_results, gsa_results, gsa_results_Dm, gsa_results_Dp)
  names(datasets) <- c("camera", "gsa_t", "gsa_Dm", "gsa_Dp")
  datasets
}

#' Return the fraction of p-values < 0.05
#' 
#' Given a Camera or GSA result, compute the fraction of p-values
#' below 0.05.
#' 
#' This function receives one of the components of the list created
#' by a call to \code{load_datasets} and returns a numeric vector
#' with the proportion of runs with a p-value below 0.05. The 
#' reason why it returns a numeric vector instead of a data frame is 
#' that this function has been created to copy and paste the output
#' into a powerpoint presentation (passing through, alas, Excel).
#' @param dataset one of the components of the list produced by 
#' calling \code{load_datasets}
#' @return a numeric vector containing, for each gene set, the 
#' fraction of p-values below 0.05.
#' @author Giovanni d'Ario
#' @export
p_signif <- function(dataset) {
  tmp <- dplyr::group_by(.data = dataset, gene_set)
  tmp <- dplyr::summarise(tmp, mean(PValue < 0.05))
  cat(tmp[[2]], sep = "\n")
}