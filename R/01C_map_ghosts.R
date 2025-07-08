map_ghosts <- function(df,
                       hydrobasin_map,
                       just_US,
                       save_maps = "figures/"){
  
  hydrobasins <- st_read(hydrobasin_map)
  
  if(just_US == TRUE){hydrobasins <- subset_hydrobasins_US(hydrobasin_map)}

  # summarize % ghosts by hydrobasin
  ghosts_by_hydro <- df %>%
    group_by(HYBAS_ID) %>%
    summarize(num_ghosts = sum(ghost, na.rm = TRUE),
              num_all = length(ghost),
              pct_ghosts = 100 * sum(ghost, na.rm = TRUE) / length(ghost),
              .groups = "drop")

  # join data to map  
  hydrobasins_ghosts <- hydrobasins %>%
    left_join(ghosts_by_hydro, by = join_by(HYBAS_ID))
  
  plot(hydrobasins_ghosts[,"pct_ghosts"])
  
  library(ggplot2)
(map_pct_ghosts <- ggplot() +
    geom_sf(data = hydrobasins_ghosts, aes(fill = pct_ghosts),
            color = "transparent") +
    scale_fill_viridis_c(option = "inferno") +
    labs(fill = "%\nmolecular\nghosts") + 
    theme_minimal() +
    theme(legend.text = element_text(size = 9)))
  
  (map_num_ghosts <-   ggplot() +
    geom_sf(data = hydrobasins_ghosts, aes(fill = num_ghosts),
            color = "transparent") +
    scale_fill_viridis_c(option = "inferno") +
    labs(fill = "No.\nmolecular\nghosts") + 
    theme_minimal() +
      theme(legend.text = element_text(size = 9)))
  
  (map_num_all <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = num_all),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "No.\nmammals") + 
      theme_minimal() +
      theme(legend.text = element_text(size = 9)))
  
  if (!is.null(save_maps)) {
    library(patchwork)
    ghost_map$num_all / ghost_map$num_ghosts / ghost_map$pct_ghosts
    ggsave(
      paste0(save_maps, "ghost_maps.png"),
      height = 5.5,
      width = 4,
      dpi = 400,
      bg = "white"
    )
  }
  
  return(list(num_ghosts = map_num_ghosts,
              num_all = map_num_all,
              pct_ghosts = map_pct_ghosts,
              data = hydrobasins_ghosts %>% as.data.frame()))
  
}