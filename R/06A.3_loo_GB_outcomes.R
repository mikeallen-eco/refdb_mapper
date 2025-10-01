loo_GB_outcomes <- function(df = loso_GB_compiled) {

# compile results for BLAST  
  df_b <- df %>%
    mutate(thresh99 = case_when(maxpident >= 99 &
                                  (maxpident - gap) < 99 &
                                  assigned_mol_name == true_mol_name ~ "correct",
                                maxpident >= 99 &
                                  (maxpident - gap) < 99 &
                                  assigned_mol_name != true_mol_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(thresh98 = case_when(maxpident >= 98 &
                           (maxpident - gap) < 98 &
                             assigned_mol_name == true_mol_name ~ "correct",
                           maxpident >= 98 &
                           (maxpident - gap) < 98 &
                             assigned_mol_name != true_mol_name ~ "incorrect",
                         TRUE ~ "abstain")) %>%
    mutate(thresh97 = case_when(maxpident >= 97 &
                                  (maxpident - gap) < 97 &
                                  assigned_mol_name == true_mol_name ~ "correct",
                                maxpident >= 97 &
                                  (maxpident - gap) < 97 &
                                  assigned_mol_name != true_mol_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(ecotag = case_when(((100 - maxpident) < gap &
                                 assigned_mol_name == true_mol_name) ~ "correct",
    ((100 - maxpident) < gap &
       assigned_mol_name != true_mol_name) ~ "incorrect",
    !((100 - maxpident) < gap) ~ "abstain"
    )) %>%
    filter(!is.na(true_mol_name)) # temporarily remove few remaining unmatched NCBI names
 
  df_b$thresh97_i <- 1*(df_b$thresh97 %in% "incorrect")
  df_b$thresh98_i <- 1*(df_b$thresh98 %in% "incorrect")
  df_b$thresh99_i <- 1 * (df_b$thresh99 %in% "incorrect")
  df_b$ecotag_i <- 1 * (df_b$ecotag %in% "incorrect")
  df_b$thresh97_c <- 1 * (df_b$thresh97 %in% "correct")
  df_b$thresh98_c <- 1 * (df_b$thresh98 %in% "correct")
  df_b$thresh99_c <- 1 * (df_b$thresh99 %in% "correct")
  df_b$ecotag_c <- 1 * (df_b$ecotag %in% "correct")
  df_b$thresh97_a <- 1 * (df_b$thresh97 %in% "abstain")
  df_b$thresh98_a <- 1 * (df_b$thresh98 %in% "abstain")
  df_b$thresh99_a <- 1 * (df_b$thresh99 %in% "abstain")
  df_b$ecotag_a <- 1 * (df_b$ecotag %in% "abstain")
 
 table(df_b$thresh98)
 table(df_b$thresh99)
 table(df_b$ecotag)
 sum(is.na(df_b$thresh98))
 
 df_b_final <- df_b %>%
   select(-c(thresh97, thresh98, thresh99, ecotag, true_ncbi_name, assigned_ncbi_name))

 return(df_b_final)
  
}