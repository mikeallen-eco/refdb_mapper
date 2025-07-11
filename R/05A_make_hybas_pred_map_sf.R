make_hybas_pred_map_sf <- function(hydrobasin_map = hybas_subset,
                                          hybas_pred = hybas_pred_i) {
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }
  
  hybas_pred_sum <- hybas_pred %>%
    group_by(HYBAS_ID) %>%
    summarize(
      mean_seqs = mean(preds_loso, na.rm = T),
      mean_ghosts = mean(preds_ghosts, na.rm = T),
      mean_all = mean(preds, na.rm = T),
      mean10 = mean(preds10, na.rm = T),
      num_above_2 = sum(preds > 2, na.rm = T),
      num10_above_2 = sum(preds10 > 2, na.rm = T),
      num_sp = length(preds),
      .groups = "drop"
    )
  
  # join hydrobasins map with predicted error data
  hydrobasins_pred <- hydrobasins %>%
    left_join(hybas_pred_sum, by = join_by(HYBAS_ID))
  
}
