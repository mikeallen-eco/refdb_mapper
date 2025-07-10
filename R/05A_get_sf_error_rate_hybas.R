get_sf_error_rate_hybas <- function(hydrobasin_map = hybas_subset,
                                    hybas_pred = hybas_pred_i) {
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }
  
  hybas_pred_sum <- hybas_pred %>%
    group_by(HYBAS_ID) %>%
    summarize(thresh99_i_mean_min1 = mean(thresh99_i, na.rm = T),
              .groups = "drop")
  
  # join hydrobasins map with predicted error data
  hydrobasins_pred <- hydrobasins %>%
    left_join(hybas_pred_sum, by = join_by(HYBAS_ID))
  
}
