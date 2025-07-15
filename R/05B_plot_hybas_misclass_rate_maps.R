plot_hybas_misclass_rate_maps <- function(hybas_pred_map_sf = hybas_pred_map_sf,
                                          to_plot = c("i", "a")){
  
  library(ggplot2)

  plot_list <- list()
  
  if("i" %in% to_plot){
  (plot_list$i <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = mean_i),
              color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "%") +
      theme_minimal())
  }
  
  if("a" %in% to_plot){
  (plot_list$a <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = mean_a),
              color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "%") +
      theme_minimal())
  }
  
  if("c" %in% to_plot){
  (plot_list$c <- hybas_pred_map_sf %>%
      ggplot() +
      geom_sf(aes(fill = mean_c),
              color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "%") +
      theme_minimal())
  }
  
  if("nnd" %in% to_plot){
    (plot_list$nnd <- hybas_pred_map_sf %>%
       ggplot() +
       geom_sf(aes(fill = mean_nnd),
               color = NA) +
       scale_fill_viridis_c(option = "inferno", direction = -1) +
       labs(fill = "MY") +
       theme_minimal())
  }
  
  return(plot_list)
  
}