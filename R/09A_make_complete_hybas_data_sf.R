# make final map sf with nnd, n_seqs, and error prediction data
make_complete_hybas_data_sf <- function(complete_hybas_data = complete_hybas_data,
                                        map = hydrobasin_map,
                                        simplify_pct = 0.01){
  
  markers <- names(complete_hybas_data)
  
  message("Simplifying map at ", simplify_pct, " ...")
  # map_simple <- map
  tictoc::tic("Done.")
  map_simple <- simplify_sf(map, p_keep = simplify_pct)
  tictoc::toc()
  
  data_base <- complete_hybas_data[[1]] %>% select(HYBAS_ID:nnd)
  
  data_by_marker <- lapply(1:length(markers), function(i) {
    data_marker <- complete_hybas_data[[i]] %>% select(n_seqs:p_c)
    colnames(data_marker) <- paste0(markers[i], "_", colnames(data_marker))
    return(data_marker)
  }) %>%
    do.call(bind_cols, .)
  
  data_final <- data_base %>%
    bind_cols(data_by_marker)
  
  map_with_data <- map_simple %>%
    left_join(data_final, 
              by = join_by(HYBAS_ID))
  
  # Fix invalid geometry in the polygon data
  message("Making sure sf has valid geometry ...")
  tictoc::tic("Done.")
  map_with_data_valid <- st_make_valid(map_with_data)
  tictoc::toc()
  
  map_with_data_valid_final <- map_with_data_valid %>%
    select(-c(MAIN_BAS:COAST))
  
  return(map_with_data_valid_final)
}
