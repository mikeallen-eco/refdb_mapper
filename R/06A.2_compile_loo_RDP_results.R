# compile GhostBLASTer & BLAST results for leave-one-sequence-out data & add info

compile_loo_RDP_results <- function(loo_path,
                                             refdb_harmonized_path = refdb_harmonized_path,
                                             verbose = F) {
  
  refdb_harmonized_df <- read.csv(refdb_harmonized_path)
  
  loo_compiled <- read_csvs(loo_path, pattern_str = "csv") %>%
    mutate(Species = ununderscore(Species)) %>%
    filter(!grepl(Species, pattern = "_x_")) %>%
    filter(!Species %in% ununderscore(ncbi_extinct)) %>%
    dplyr::rename(true_ncbi_name = target_sp) %>%
    mutate(true_ncbi_name = ununderscore(true_ncbi_name)) %>%
    select(-acc) %>%
    filter(!grepl(qseqid, pattern = "_x_")) %>%
    filter(!true_ncbi_name %in% ununderscore(ncbi_extinct)) %>%
    left_join(refdb_harmonized_df %>% select(true_ncbi_name = full_sci_name,
                                          true_mol_name = BB_Accepted),
              by = join_by(true_ncbi_name)) %>%
    mutate(assigned_ncbi_name = Species) %>%
    left_join(refdb_harmonized_df %>% select(assigned_ncbi_name = full_sci_name,
                                          assigned_mol_name = BB_Accepted),
              by = join_by(assigned_ncbi_name)) %>%
    select(qseqid, true_mol_name, assigned_mol_name, 
           spboot, 
           true_ncbi_name, assigned_ncbi_name)
  
  return(loo_compiled)
}