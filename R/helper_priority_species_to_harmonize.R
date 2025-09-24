# find widespread species that still need manual taxonomic matching (to prioritize)
priority_species_to_harmonize <- function(harmonized_df, 
                                          hydrobasin_df = hydrobasin_refdb_info,
                                          show = 5) {
  to_go <- harmonized_df[is.na(harmonized_df$BB_Accepted),]$uid %>%
    as.data.frame() %>% dplyr::rename(sp = 1)
  
  sort(table(hydrobasin_df$mol_name), decreasing = T) %>% 
    as.data.frame() %>%
    dplyr::rename(sp = 1) %>%
    mutate(to_go = case_when(sp %in% to_go$sp ~ 1,
                             TRUE ~ 0)) %>% arrange(desc(to_go), desc(Freq)) %>% head(n = show)
}
