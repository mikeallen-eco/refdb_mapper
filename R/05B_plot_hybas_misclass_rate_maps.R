plot_hybas_misclass_rate_maps <- function(hybas_pred_map_sf = hybas_pred_map_sf,
                                          save_to = "figures/",
                                          suffix = "",
                                          h = 6, w = 8, res = 400){
  
  library(ggplot2)
  library(patchwork)
  
  (current <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = mean_all),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno", limits = c(0,8.5)) +
      labs(fill = "Mean\npredicted\nmisclass.\nrate (%)",
           title = "Current reference database\n(51% ghosts, med. 1.5 seqs / sp)") +
      theme_minimal())
  
  (future <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = mean10),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno", limits = c(0,8)) +
      labs(fill = "Mean\npredicted\nmisclass.\nrate (%)",
           title = "Future reference database\n(≥10 seqs / sp)") +
      theme_minimal())
  
  (num_current <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = log10(num_above_2)),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno", limits = c(-0.35,2.3),
                           breaks = log10(c(1,3,10,30,100)), labels = c(1,3,10,30,100)) +
      labs(fill = "No. sp. with\npredicted\nmisclass.\nrate >2%",
           title = "Current reference database\n(51% ghosts, med. 1.5 seqs / sp)") +
      theme_minimal())
  
  (num_future <- hybas_pred_map_sf %>%
      mutate(num10_above_2 = num10_above_2 + 0.5) %>%
      ggplot() +
      geom_sf(aes(fill = log10(num10_above_2)),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno", limits = c(-0.35,2.3),
                           breaks = log10(c(1,3,10,30,100)), labels = c(1,3,10,30,100)) +
      labs(fill = "No. sp. with\npredicted\nmisclass.\nrate >2%",
           title = "Future reference database\n(≥10 seqs / sp)") +
      theme_minimal())
  
  # make composite plot
  composite_plot <- (current | future) / (num_current | num_future)
  
  if(!is.null(save_to)){
    ggsave(
      paste0(save_to, "current_future_misclass_map_", suffix, ".png"),
      plot = composite_plot,
      height = h,
      width = w,
      dpi = res,
      bg = "white"
    )
  }
  
  return(list(current = current,
              future = future,
              num_current = num_current,
              num_future = num_future))
  
}