plot_predicted_loso_lopso_error <- function(preds = preds_loso_lospo) {
  
  library(ggplot2)

  LOSO_plot <- plot_predicted_LOSO_error(preds$preds_gb_loso)
  LOSpO_plot <- plot_predicted_LOSpO_error(preds$preds_gb_lospo)

  return(list(loso = LOSO_plot,
              lospo = LOSpO_plot))
  
}