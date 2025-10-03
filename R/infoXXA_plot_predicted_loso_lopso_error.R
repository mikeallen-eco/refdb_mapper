library(ggplot2)

plot_predicted_loso_lopso_error <- function(preds = fits,
                                            markers = markers,
                                            rubrics = rubrics) {
  
  rubrics_list <- lapply(rubrics, function (j){
    
  preds_rubric <- preds[[j]]
  
  marker_list <- lapply(1:length(markers), function (i){
  
  LOSO_plots <- plot_predicted_LOSO_error(preds_loso = preds_rubric[[i]]$loso)
  LOSpO_plots <- plot_predicted_LOSpO_error(preds_lospo = preds_rubric[[i]]$lospo)

  return(list(loso = LOSO_plots,
              lospo = LOSpO_plots))
  })
  
  names(marker_list) <- markers
  return(marker_list)
  })
  
  names(rubrics_list) <- rubrics
  
  return(rubrics_list)
}

