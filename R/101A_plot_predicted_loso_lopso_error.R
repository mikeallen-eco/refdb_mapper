library(ggplot2)

plot_predicted_loso_lopso_error <- function(preds = fits,
                                            markers = markers) {
  
  marker_list <- lapply(1:length(markers), function (i){

  LOSO_plots <- plot_predicted_LOSO_error(preds$loso[[i]])
  LOSpO_plots <- plot_predicted_LOSpO_error(preds$lospo[[i]])

  return(list(loso = LOSO_plots,
              lospo = LOSpO_plots))
  })
  
  names(marker_list) <- markers
  
  return(marker_list)
}