summarize_by_hydrobasin <- function(df = hydrobasin_ghosts){

# summarize % ghosts by hydrobasin
ghosts_by_hydro <- df %>%
  group_by(HYBAS_ID) %>%
  summarize(num_ghosts = sum(ghost, na.rm = TRUE),
            num_all = length(ghost),
            pct_ghosts = 100 * sum(ghost, na.rm = TRUE) / length(ghost),
            .groups = "drop")

medians <- df %>%
  filter(n_seqs > 0) %>%
  group_by(HYBAS_ID) %>%
  summarize(med_seqs = median(n_seqs, na.rm = TRUE),
            .groups = "drop")

ghosts_by_hydro_w_medians <- ghosts_by_hydro %>%
  left_join(medians %>% select(HYBAS_ID, med_seqs), 
            by = join_by(HYBAS_ID))

return(ghosts_by_hydro_w_medians)

}