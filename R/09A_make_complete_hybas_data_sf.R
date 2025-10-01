# make final map sf with nnd, n_seqs, and error prediction data
make_complete_hybas_data_sf <- function(complete_hybas_preds = list(hybas_preds_blast97,
                                                                    hybas_preds_blast98,
                                                                    hybas_preds_blast99),
                                        map = hydrobasin_map,
                                        simplify_pct = 0.01){
  
  markers <- names(complete_hybas_preds[[1]])
  
  message("Simplifying map at ", simplify_pct, " ...")
  # map_simple <- map
  tictoc::tic("Done.")
  map_simple <- simplify_sf(map, p_keep = simplify_pct)
  tictoc::toc()
  
  markers <- names(complete_hybas_preds[[1]])
  
  data_base <- complete_hybas_preds[[1]][[1]] %>% select(HYBAS_ID:phyl_name)
  
  rubric_list <- list()
  for(x in 1:length(complete_hybas_preds)){
    
  data_by_rubric <- lapply(1:length(markers), function(i) {
    data_marker <- complete_hybas_preds[[x]][[i]] %>% select(contains(markers[i]))
    # colnames(data_marker) <- paste0(markers[i], "_", colnames(data_marker))
    return(data_marker)
  }) %>%
    do.call(bind_cols, .)
  
  rubric_list[[x]] <- data_by_rubric
  
  }
  
  rubrics_df <- rubric_list %>%
    do.call(bind_cols, .)
  
  data_final <- data_base %>%
    bind_cols(rubrics_df)
  
  map_with_data <- map_simple %>%
    select(-c(MAIN_BAS:COAST)) %>%
    left_join(data_final, 
              by = join_by(HYBAS_ID))
  
  # Fix invalid geometry in the polygon data
  message("Making sure sf has valid geometry ...")
  tictoc::tic("Done.")
  map_with_data_valid_final <- st_make_valid(map_with_data)
  tictoc::toc()
  
  return(map_with_data_valid_final)
}
