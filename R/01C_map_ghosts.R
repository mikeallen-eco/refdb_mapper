map_ghosts <- function(df,
                       hydrobasin_map,
                       save_maps = "figures/"){
  
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }  

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
    map_num_all / map_num_ghosts / map_pct_ghosts
    ggsave(
      paste0(save_maps, "ghost_maps_NA.png"),
      height = 7.5,
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