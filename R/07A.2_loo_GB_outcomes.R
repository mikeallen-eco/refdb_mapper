loo_GB_outcomes <- function(df) {

# compile results for BLAST  
  df_b <- df %>%
    mutate(thresh99 = case_when(max_pident_nei_no_ecto >= 99 &
                                  (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 99 &
                                  top_match_nei_no_ecto == true_ncbi_name ~ "correct",
                                max_pident_nei_no_ecto >= 99 &
                                  (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 99 &
                                  top_match_nei_no_ecto != true_ncbi_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(thresh98 = case_when(max_pident_nei_no_ecto >= 98 &
                           (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 98 &
                           top_match_nei_no_ecto == true_ncbi_name ~ "correct",
                         max_pident_nei_no_ecto >= 98 &
                           (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 98 &
                           top_match_nei_no_ecto != true_ncbi_name ~ "incorrect",
                         TRUE ~ "abstain")) %>%
    mutate(ecotag = case_when(((100 - max_pident_nei_no_ecto) < max_pident_gap_nei_no_ecto &
                                 top_match_nei_no_ecto == true_ncbi_name) ~ "correct",
    ((100 - max_pident_nei_no_ecto) < max_pident_gap_nei_no_ecto &
       top_match_nei_no_ecto != true_ncbi_name) ~ "incorrect",
    !((100 - max_pident_nei_no_ecto) < max_pident_gap_nei_no_ecto) ~ "abstain"
    )) %>%
    mutate(conf90 = case_when((score_nei_no_ecto >= 90 &
                                 top_match_nei_no_ecto == true_ncbi_name) ~ "correct",
                              (score_nei_no_ecto >= 90 &
                                 top_match_nei_no_ecto != true_ncbi_name) ~ "incorrect",
                              score_nei_no_ecto < 90 ~ "abstain")) %>%
    mutate(method = "BLAST")
  
  # compile results for GhostBLASTer
  df_gb <- df %>%
    mutate(thresh99 = case_when(max_pident_nei_ecto >= 99 &
                                  (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 99 &
                                  top_match_nei_ecto == true_ncbi_name ~ "correct",
                                max_pident_nei_ecto >= 99 &
                                  (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 99 &
                                  top_match_nei_ecto != true_ncbi_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(thresh98 = case_when(max_pident_nei_ecto >= 98 &
                                  (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 98 &
                                  top_match_nei_ecto == true_ncbi_name ~ "correct",
                                max_pident_nei_ecto >= 98 &
                                  (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 98 &
                                  top_match_nei_ecto != true_ncbi_name ~ "incorrect",
                                TRUE ~ "abstain")) %>%
    mutate(ecotag = case_when(((100 - max_pident_nei_ecto) < max_pident_gap_nei_ecto &
                                 top_match_nei_ecto == true_ncbi_name) ~ "correct",
                              ((100 - max_pident_nei_ecto) < max_pident_gap_nei_ecto &
                                 top_match_nei_ecto != true_ncbi_name) ~ "incorrect",
                              !((100 - max_pident_nei_ecto) < max_pident_gap_nei_ecto) ~ "abstain"
    )) %>%
    mutate(conf90 = case_when((score_nei_ecto >= 90 &
                                 top_match_nei_ecto == true_ncbi_name) ~ "correct",
                              (score_nei_ecto >= 90 &
                                 top_match_nei_ecto != true_ncbi_name) ~ "incorrect",
                              score_nei_ecto < 90 ~ "abstain")) %>%
    mutate(method = "GhostBLASTer")
  
  # select & format final df form
 df_both <- df_b %>%
   bind_rows(df_gb)
 
 df_both$thresh98_i <- 1*(df_both$thresh98 %in% "incorrect")
 df_both$thresh99_i <- 1*(df_both$thresh99 %in% "incorrect")
 df_both$ecotag_i <- 1*(df_both$ecotag %in% "incorrect")
 df_both$conf90_i <- 1*(df_both$conf90 %in% "incorrect")
 df_both$thresh98_c <- 1*(df_both$thresh98 %in% "correct")
 df_both$thresh99_c <- 1*(df_both$thresh99 %in% "correct")
 df_both$ecotag_c <- 1*(df_both$ecotag %in% "correct")
 df_both$conf90_c <- 1*(df_both$conf90 %in% "correct")
 df_both$thresh98_a <- 1*(df_both$thresh98 %in% "abstain")
 df_both$thresh99_a <- 1*(df_both$thresh99 %in% "abstain")
 df_both$ecotag_a <- 1*(df_both$ecotag %in% "abstain")
 df_both$conf90_a <- 1*(df_both$conf90 %in% "abstain")
 
 df_both <- df_both %>%
   select(qseqid, method, starts_with("top_match"), 
          starts_with("score"), starts_with("max"), 
          true_ncbi_name:conf90_a) %>%
   dplyr::rename(n_seqs = n)
 
 # table(df_both$method, df_both$thresh98)
 # table(df_both$method, df_both$thresh99)
 # sum(is.na(df_both$thresh99))
return(df_both)
  
}