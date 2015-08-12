#' Plot all the main results from the simulations
#' 
#' Go through the folders produced by running the simulations and
#' create the plots for each of them.
#' 
#' This function checks if there are folders starting with 'fc'
#' and reads the content of their subdirectories, creating a plot
#' for each of the comparisons stored therein. It creates then a
#' plot folder and, for each subfolder, generates a subfolder with
#' the relative plots
#' @return No return value but, as a side effect, plots are created.
#' @author Giovanni d'Ario
#' @export
plot_results <- function() {
  require(ggplot2)
  if (!file.exists('plots'))
    dir.create('plots')
  
  plot_results <- ggplot() + aes(x = PValue, fill = gene_set > 5) +
    geom_histogram() +  
    scale_fill_discrete(labels = c("Yes", "No")) + 
    guides(fill = guide_legend(title = "Enriched")) +
    facet_wrap(~ gene_set, nrow = 2)
  
  dirs <- dir(pattern = "^fc_")
  
  for(d in dirs) {
    subdirs <- dir(d, full.names = TRUE)
    for(subd in subdirs) {
      
      dir.create(file.path("plots", subd), recursive = TRUE)
      plotdir <- file.path("plots", subd)
      
      load(file.path(subd, "camera_results.RData"))
      
      p1 <- plot_results %+% camera_results
      ggsave(p1, filename = file.path(plotdir, "camera_results.png"), 
        width = 10, height = 8)
      
      load(file.path(subd, "gsa_results.RData"))
      
      p2 <- plot_results %+% gsa_results
      ggsave(p2, filename = file.path(plotdir, "gsa_results_t_test.png"), 
        width = 10, height = 8)
      
      load(file.path(subd, "gsa_results_KSw_Dm.RData"))
      
      p3 <- plot_results %+% gsa_results_Dm
      ggsave(p3, filename = file.path(plotdir, "gsa_results_KSw_Dm.png"), 
        width = 10, height = 8)
      
      load(file.path(subd, "gsa_results_KSw_Dp.RData"))
      
      p4 <- plot_results %+% gsa_results_Dp
      ggsave(p4, filename = file.path(plotdir, "gsa_results_KSw_Dp.png"), 
        width = 10, height = 8)
    }
  }
}