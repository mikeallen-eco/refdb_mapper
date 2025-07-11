plot_predicted_loso_lopso_error <- function(preds = preds_loso_lospo,
                                            save_to = "figures/",
                                            suffix = "",
                                            h = 4,
                                            w = 8,
                                            res = 400) {
  library(ggplot2)
  library(patchwork)
  LOSO_plot <- plot_predicted_LOSO_error(preds$preds_gb_loso)
  LOSpO_plot <- plot_predicted_LOSpO_error(preds$preds_gb_lospo)
  combined_plot_i <- LOSO_plot$i | LOSpO_plot$i
  combined_plot_a <- LOSO_plot$a | LOSpO_plot$a
  combined_plot_c <- LOSO_plot$c

  if (!is.null(save_to)) {
    ggsave(
      paste0(
        save_to,
        "LOSO_LOSpO_predicted_misclassified_",
        suffix,
        ".png"
      ),
      plot = combined_plot_i,
      height = h,
      width = w,
      dpi = res,
      bg = "white"
    )
    
    ggsave(
      paste0(
        save_to,
        "LOSO_LOSpO_predicted_unclassified_",
        suffix,
        ".png"
      ),
      plot = combined_plot_a,
      height = h,
      width = w,
      dpi = res,
      bg = "white"
    )
    
    ggsave(
      paste0(
        save_to,
        "LOSO_LOSpO_predicted_%correct_",
        suffix,
        ".png"
      ),
      plot = combined_plot_a,
      height = h,
      width = w/2,
      dpi = res,
      bg = "white"
    )
  }
  
  return(list(i = combined_plot_i,
              a = combined_plot_a,
              c = combined_plot_c))
  
}