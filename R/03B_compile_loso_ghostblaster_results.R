# compile GhostBLASTer & BLAST results for leave-one-sequence-out data & add info

compile_loso_ghostblaster_results <- function(loso_path = paste0(out_path, "loso"),
                                       loso_refdb_ndd_df = loso_refdb_ndd,
                                       verbose = FALSE) {
  
  loso_compiled <- read_csvs(loso_path, pattern_str = "csv") %>%
    mutate(local = 1) %>% 
    ghostsum(., marker = "Vences_16S", verbose = verbose) %>%
    mutate(tmp = qseqid) %>%
    separate(tmp, into = c("g", "s", "acc"), sep = "_") %>%
    mutate(true_ncbi_name = paste0(g, " ", s)) %>%
    select(-g, -s, acc) %>%
    filter(true_ncbi_name %in% loso_refdb_ndd$ncbi_name) %>%
    left_join(loso_refdb_ndd %>% dplyr::rename(true_ncbi_name = ncbi_name),
              by = join_by(true_ncbi_name)) %>%
    select(-phyl_name) %>%
    filter(!grepl(top_match_loc_ecto, pattern = "Skip"))
  
  return(loso_compiled)
}