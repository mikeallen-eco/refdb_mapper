map_ghosts <- function(df = seq_info_summarized_by_hydrobasin,
                       hydrobasin_map,
                       to_plot = c("num_ghosts", "pct_ghosts", 
                                   "num_all")){
  
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }  

  # join data to map  
  hydrobasins_ghosts <- hydrobasins %>%
    left_join(df, by = join_by(HYBAS_ID))
  
  plot_list <- list()
  
  library(ggplot2)
  if ("pct_ghosts" %in% to_plot) {
    plot_list$pct_ghosts <- ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = pct_ghosts), color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "%\nmolecular\nghosts") +
      theme_minimal() +
      theme(legend.text = element_text(size = 9))
  }

  if ("num_ghosts" %in% to_plot) {
    plot_list$num_ghosts <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = num_ghosts), color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "No.\nmolecular\nghosts") +
      theme_minimal() +
      theme(legend.text = element_text(size = 9))
  }
  
  if ("num_all" %in% to_plot) {
    plot_list$num_all <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = num_all), color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "No.\nmammals") +
      theme_minimal() +
      theme(legend.text = element_text(size = 9))
  }
  
  if ("med_seqs" %in% to_plot) {
  plot_list$med_seqs <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = med_seqs),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno", limits = c(0,7)) +
      labs(fill = "Median\nnum. seqs.\n(non-ghosts)") + 
      theme_minimal() +
      theme(legend.text = element_text(size = 9))
  }
  
  return(plot_list)
  
}