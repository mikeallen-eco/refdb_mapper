# compile GhostBLASTer & BLAST results for leave-one-sequence-out data & add info

compile_loo_ghostblaster_results <- function(loo_path,
                                             refdb_harmonized_path = refdb_harmonized_path,
                                             verbose = F) {
  
  refdb_harmonized_df <- read.csv(refdb_harmonized_path)
  
  loo_compiled <- read_csvs(loo_path, pattern_str = "csv") %>%
    select(-multi_class) %>%
    patch_old_gb_version() %>%
    mutate(local = 1) %>%
    filter(!grepl(seq_species, pattern = "_x_")) %>%
    filter(!seq_species %in% ununderscore(ncbi_extinct)) %>%
    ghostsum(., marker = "average", verbose = verbose) %>%
    mutate(tmp = qseqid) %>%
    separate(tmp, into = c("g", "s", "acc"), sep = "_") %>%
    mutate(true_ncbi_name = paste0(g, " ", s)) %>%
    select(-g, -s, acc) %>%
    filter(!grepl(qseqid, pattern = "_x_")) %>%
    filter(!true_ncbi_name %in% ununderscore(ncbi_extinct)) %>%
    mutate(maxpident = case_when(grepl(blastg_local, pattern = "Skip") ~ 75,
                                           TRUE ~ maxpident_blast_all),
           gap = case_when(grepl(blastg_local, pattern = "Skip") ~ 0,
                                               TRUE ~ gap_blast_all),
           blast_all = case_when(grepl(blastg_local, pattern = "Skip") ~ "skipped",
                                          TRUE ~ blast_all)) %>%
    left_join(refdb_harmonized_df %>% select(true_ncbi_name = full_sci_name,
                                          true_mol_name = BB_Accepted),
              by = join_by(true_ncbi_name)) %>%
    mutate(assigned_ncbi_name = blast_all) %>%
    left_join(refdb_harmonized_df %>% select(assigned_ncbi_name = full_sci_name,
                                          assigned_mol_name = BB_Accepted),
              by = join_by(assigned_ncbi_name)) %>%
    mutate(assigned_mol_name = case_when(grepl(assigned_ncbi_name, pattern = "skip") ~ "skipped",
                                         TRUE ~ assigned_mol_name)) %>%
    select(qseqid, true_mol_name, assigned_mol_name, 
           maxpident, gap, 
           true_ncbi_name, assigned_ncbi_name)
  
  return(loo_compiled)
}