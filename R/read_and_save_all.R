#' Read the files created by Camera and GSA
#' 
#' Read the files created by Camera and GSA, and save them as 
#' data frames where the information is (mostly) the same for
#' the output coming from the two methods.
#' 
#' This function identifies all the folders starting with 'fc_'
#' and reads the files produced by the simulation storing the 
#' relevant information in data frames that are then saved in
#' the same folder.
#' @return no return value but, as a side effect, a few RData 
#' files are created
#' @author Giovanni d'Ario
#' @export
read_and_save_all <- function() {
  
  results_dirs <- dir(pattern = "^fc_")
  
  for(resdir in results_dirs) {
    results_subdirs <- dir(resdir, full.names = TRUE)
    for(ressubdir in results_subdirs) {
      message(paste("processing", ressubdir))
      
      ## Camera results
      tmp <- read_camera_results(resultsFile = file.path(ressubdir, 
        "results_simulation.RData"), outFile = "camera_results.RData")
      
      ## Save all the results from gsa
      tmp <- save_gsa_results(dirName = ressubdir, 
        outFile = "gsa_results_all_tests.RData")
      
      ## Extract the KSw Dp results from gsa_results
      gsa_results_Dp <- read_gsa_results(file.path(ressubdir, 
        "gsa_results_all_tests.RData"), stat = "fc_KSw_Dp")
      
      ## Save the Dp results
      save(gsa_results_Dp, file = file.path(ressubdir, "gsa_results_KSw_Dp.RData"))
      
      ## Extract the KSw Dm results from gsa_results
      gsa_results_Dm <- read_gsa_results(file.path(ressubdir, 
        "gsa_results_all_tests.RData"), stat = "fc_KSw_Dm")
      
      ## Save the Dm results
      save(gsa_results_Dm, file = file.path(ressubdir, "gsa_results_KSw_Dm.RData"))
    }
  }
}