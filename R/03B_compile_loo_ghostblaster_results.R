# compile GhostBLASTer & BLAST results for leave-one-sequence-out data & add info

compile_loo_ghostblaster_results <- function(loo_path,
                                       loo_refdb_ndd_df,
                                       include_skips = TRUE,
                                       verbose = FALSE) {
  
  loo_compiled <- read_csvs(loo_path, pattern_str = "csv") %>%
    mutate(local = 1) %>% 
    ghostsum(., marker = "Vences_16S", verbose = verbose) %>%
    mutate(tmp = qseqid) %>%
    separate(tmp, into = c("g", "s", "acc"), sep = "_") %>%
    mutate(true_ncbi_name = paste0(g, " ", s)) %>%
    select(-g, -s, acc) %>%
    filter(true_ncbi_name %in% loo_refdb_ndd_df$ncbi_name) %>%
    left_join(loo_refdb_ndd_df %>% dplyr::rename(true_ncbi_name = ncbi_name),
              by = join_by(true_ncbi_name)) %>%
    select(-phyl_name) %>%
    mutate(max_pident_nei_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 75,
                                               TRUE ~ max_pident_nei_ecto),
           max_pident_gap_nei_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 0,
                                                   TRUE ~ max_pident_gap_nei_ecto),
           max_pident_nei_no_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 75,
                                                  TRUE ~ max_pident_nei_no_ecto),
           max_pident_gap_nei_no_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 0,
                                                      TRUE ~ max_pident_gap_nei_no_ecto),
           top_match_nei_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ "skipped",
                                              TRUE ~ top_match_nei_ecto),
           top_match_nei_no_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ "skipped",
                                                 TRUE ~ top_match_nei_no_ecto),
           score_nei_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 0,
                                          TRUE ~ score_nei_ecto),
           score_nei_no_ecto = case_when(grepl(top_match_loc_ecto, pattern = "Skip") ~ 0,
                                             TRUE ~ score_nei_no_ecto))
  
  if(include_skips == FALSE){
    loo_compiled <- loso_compiled %>%
    filter(!grepl(top_match_loc_ecto, pattern = "Skip"))
  }
  
  return(loo_compiled)
}