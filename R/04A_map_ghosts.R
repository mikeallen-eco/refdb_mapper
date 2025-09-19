library(ggplot2)

map_ghosts <- function(df = hydrobasin_ref_nnd_info_sum,
                       hydrobasin_map,
                       markers = markers,
                       to_plot = c("num_ghosts", "pct_ghosts", 
                                   "total_species")){
  
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }  
  
  # loop through markers
  final_list <- lapply(1:length(markers), function(i){
  
  plot_list <- list()
  
  # subset to 1 marker
  df2 <- df %>%
    select(HYBAS_ID, contains(markers[i]), total_species, nnd_med) %>%
    rename_with(~"pct_ghosts", contains("pct_ghosts")) %>%
    rename_with(~"num_ghosts", contains("num_ghosts")) %>%
    rename_with(~"med_seqs", contains("med_seqs"))
    
  # join data to map  
  hydrobasins_ghosts <- hydrobasins %>% 
    left_join(df2, by = join_by(HYBAS_ID))
                                
  if ("pct_ghosts" %in% to_plot) {
    plot_list$pct_ghosts <- ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = pct_ghosts), color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "%\nmolecular\nghosts") +
      theme_minimal() +
      theme(legend.text = element_text(size = 9)) +
      ggtitle(markers[i])
  }

  if ("num_ghosts" %in% to_plot) {
    plot_list$num_ghosts <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = num_ghosts), color = NA) +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "No.\nmolecular\nghosts") +
      theme_minimal() +
      theme(legend.text = element_text(size = 9)) +
      ggtitle(markers[i])
  }
  
  if ("total_species" %in% to_plot) {
    plot_list$total_species <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = total_species), color = NA) +
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
      theme(legend.text = element_text(size = 9)) +
    ggtitle(markers[i])
  }
  
  if ("nnd_med" %in% to_plot) {
    plot_list$nnd_med <-   ggplot() +
      geom_sf(data = hydrobasins_ghosts, aes(fill = nnd_med),
              color = "transparent") +
      scale_fill_viridis_c(option = "inferno") +
      labs(fill = "Median\nNND\nw/in basin") + 
      theme_minimal() +
      theme(legend.text = element_text(size = 9))
  }

  message(markers[i])
  
  return(plot_list)
  })

  names(final_list) <- markers
  return(final_list)
  }