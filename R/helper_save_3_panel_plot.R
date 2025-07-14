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