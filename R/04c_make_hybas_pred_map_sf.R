make_hybas_pred_map_sf <- function(hydrobasin_map = hybas_map,
                                   hybas_pred = hybas_error_data,
                                   ghost_info = seq_info_summarized_by_hydrobasin) {
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }
  
  hybas_pred_sum <- hybas_pred %>%
    group_by(HYBAS_ID) %>%
    summarize(
      mean_nnd = mean(nnd, na.rm = T),
      mean_i = mean(preds_i, na.rm = T),
      mean_a = mean(preds_a, na.rm = T),
      mean_c = mean(preds_c, na.rm = T),
      mean10_i = mean(preds10_i, na.rm = T),
      mean10_a = mean(preds10_a, na.rm = T),
      mean10_c = mean(preds10_c, na.rm = T),
      num_above_2_i = sum(preds_i > 2, na.rm = T),
      num10_above_2_i = sum(preds10_i > 2, na.rm = T),
      num_sp = length(preds_i),
      .groups = "drop"
    )
  
  # join hydrobasins map with predicted error data
  hydrobasins_pred <- hydrobasins %>%
    left_join(hybas_pred_sum, by = join_by(HYBAS_ID)) %>%
    left_join(ghost_info, 
              by = join_by(HYBAS_ID))
  
}
