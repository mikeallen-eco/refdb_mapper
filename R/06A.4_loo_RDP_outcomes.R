loo_RDP_outcomes <- function(df = loso_RDP_compiled) {

# compile results for BLAST  
  df_rdp <- df %>%
    mutate(rdp70 = case_when(spboot >= 70 &
                               assigned_mol_name == true_mol_name ~ "correct",
                             spboot >= 70 &
                               assigned_mol_name != true_mol_name ~ "incorrect",
                             TRUE ~ "abstain")) %>%
    mutate(rdp80 = case_when(spboot >= 80 &
                                  assigned_mol_name == true_mol_name ~ "correct",
                             spboot >= 80 &
                                  assigned_mol_name != true_mol_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(rdp90 = case_when(spboot >= 90 &
                                  assigned_mol_name == true_mol_name ~ "correct",
                                spboot >= 90 &
                                  assigned_mol_name != true_mol_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(rdp95 = case_when(spboot >= 95 &
                                  assigned_mol_name == true_mol_name ~ "correct",
                                spboot >= 95 &
                                  assigned_mol_name != true_mol_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    filter(!is.na(true_mol_name)) # temporarily remove few remaining unmatched NCBI names
 
  df_rdp$rdp70_i <- 1*(df_rdp$rdp70 %in% "incorrect")
  df_rdp$rdp80_i <- 1*(df_rdp$rdp80 %in% "incorrect")
  df_rdp$rdp90_i <- 1 * (df_rdp$rdp90 %in% "incorrect")
  df_rdp$rdp95_i <- 1 * (df_rdp$rdp95 %in% "incorrect")
  df_rdp$rdp70_c <- 1*(df_rdp$rdp70 %in% "correct")
  df_rdp$rdp80_c <- 1*(df_rdp$rdp80 %in% "correct")
  df_rdp$rdp90_c <- 1 * (df_rdp$rdp90 %in% "correct")
  df_rdp$rdp95_c <- 1 * (df_rdp$rdp95 %in% "correct")
  df_rdp$rdp70_a <- 1*(df_rdp$rdp70 %in% "abstain")
  df_rdp$rdp80_a <- 1*(df_rdp$rdp80 %in% "abstain")
  df_rdp$rdp90_a <- 1 * (df_rdp$rdp90 %in% "abstain")
  df_rdp$rdp95_a <- 1 * (df_rdp$rdp95 %in% "abstain")
 
 df_rdp_final <- df_rdp %>%
   select(-c(rdp70, rdp80, rdp90, rdp95, true_ncbi_name, assigned_ncbi_name))

 return(df_rdp_final)
  
}