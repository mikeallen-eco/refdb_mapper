error_rate_scatterplots <- function(hybas_pred_map_sf = hybas_pred_map_sf){ 
  (i <- hybas_pred_map_sf %>%
     as.data.frame() %>%
     ggplot() +
     geom_point(aes(x = pct_ghosts, y = mean_i, color = mean_nnd), alpha = 0.7) +
     geom_point(aes(x = pct_ghosts, y = mean_i), data = hybas_pred_map_sf %>% as.data.frame() %>% filter(num_all >150),
                shape = "+", color = "black", size = 1) +
     scale_color_viridis_c() +
     theme_minimal() +
     labs(x = "% molecular ghosts", 
          y = "Mean predicted rate (%)",
          color = "Mean\nNND\n(MY)",
          title = "Misclassification"))
  
  (a <- hybas_pred_map_sf %>%
      as.data.frame() %>%
      ggplot() +
      geom_point(aes(x = pct_ghosts, y = mean_a, color = mean_nnd), alpha = 0.7) +
      geom_point(aes(x = pct_ghosts, y = mean_a), data = hybas_pred_map_sf %>% as.data.frame() %>% filter(num_all >150),
                 shape = "+", color = "black", size = 1) +
      scale_color_viridis_c() +
      theme_minimal() +
      labs(x = "% molecular ghosts", 
           y = "Mean predicted rate (%)",
           color = "Mean\nNND\n(MY)",
           title = "Unclassified"))
  
  (c <- hybas_pred_map_sf %>%
      as.data.frame() %>%
      ggplot() +
      geom_point(aes(x = pct_ghosts, y = mean_c, color = mean_nnd), alpha = 0.7) +
      geom_point(aes(x = pct_ghosts, y = mean_c), data = hybas_pred_map_sf %>% as.data.frame() %>% filter(num_all >150),
                 shape = "+", color = "black", size = 1) +
      scale_color_viridis_c() +
      theme_minimal() +
      labs(x = "% molecular ghosts", 
           y = "Mean predicted rate (%)",
           color = "Mean\nNND\n(MY)",
           title = "Correct Classification"))
  
  return(list(i = i, a = a, c = c))
}