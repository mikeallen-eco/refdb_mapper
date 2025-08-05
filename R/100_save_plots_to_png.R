# functions to save various plots to png
library(ggplot2)
library(patchwork)

save_num_mammals <- function(plot = ghost_plots$num_all) {
  ggsave(
    "figures/num_mammals.png",
    bg = "white",
    height = 6,
    width = 9,
    plot = plot,
    dpi = 400
  )
}

save_mean_nnd <- function(plot = hybas_mean_nnd_map$nnd) {
  ggsave(
    "figures/mean_nnd.png",
    bg = "white",
    height = 6,
    width = 9,
    plot = plot,
    dpi = 400
  )
}

save_num_pct_molecular_ghosts <- function(plot = ghost_plots$num_ghosts /
                                            ghost_plots$pct_ghosts) {
  ggsave(
    "figures/num_pct_molecular_ghosts.png",
    bg = "white",
    height = 12,
    width = 9,
    plot = plot,
    dpi = 400
  )
}

save_median_num_seqs_non_ghosts <- function(plot = med_seqs$med_seqs) {
  ggsave(
    "figures/median_num_seqs_non_ghosts.png",
    bg = "white",
    height = 6,
    width = 9,
    plot = plot,
    dpi = 400
  )
}

save_3_panel_plot <- function(plot_list,
                              save_to = "figures/plot.png",
                              h = 4, w = 8, res = 400){
  library(ggplot2)
  library(patchwork)
  
  (plot_list[[1]]) / (plot_list[[2]] | plot_list[[3]]) + plot_layout(heights = c(1,0.5, 0.5))
  
  ggsave(
    save_to,
    height = h,
    width = w,
    dpi = 400,
    bg = "white"
  )
}

save_predicted_pct_misclassified_3panel <- function(plot_list = list(
  hybas_misclass_rate_maps$i,
  predicted_loso_lospo_error_plots$loso$i,
  predicted_loso_lospo_error_plots$lospo$i
)) {
  save_3_panel_plot(
    plot_list = plot_list,
    save_to = "figures/predicted_pct_misclassified_3panel.png",
    h = 10,
    w = 10,
    res = 400
  )
}

save_predicted_pct_unclassified_3panel <- function(plot_list = list(
  hybas_misclass_rate_maps$a,
  predicted_loso_lospo_error_plots$loso$a,
  predicted_loso_lospo_error_plots$lospo$a
)) {
  save_3_panel_plot(
    plot_list = plot_list,
    save_to = "figures/predicted_pct_unclassified_3panel.png",
    h = 10,
    w = 10,
    res = 400
  )
}

