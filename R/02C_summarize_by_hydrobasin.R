summarize_by_hydrobasin <- function(df = hydrobasin_ghosts_all){

  names(df)
# summarize % ghosts by hydrobasin
  ghosts_by_hydro <- df %>%
    group_by(HYBAS_ID) %>%
    summarise(across(
      # just drop the taxonomic text columns
      .cols = -c(sciname, order, family),
      .fns = list(
        num_ghosts = ~ sum(.x > 0, na.rm = TRUE),
        num_all    = ~ sum(!is.na(.x)),
        pct_ghosts = ~ 100 * sum(.x > 0, na.rm = TRUE) / sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ), .groups = "drop")

  medians <- df %>%
    group_by(HYBAS_ID) %>%
    summarise(across(
      .cols = -c(sciname, order, family),
      .fns  = ~ median(.x[.x > 0], na.rm = TRUE),  # median of positive seq counts
      .names = "{.col}_med_seqs"
    ), .groups = "drop")
  
  ghosts_by_hydro_w_medians <- ghosts_by_hydro %>%
    left_join(medians, by = join_by(HYBAS_ID))
  
return(ghosts_by_hydro_w_medians)

}