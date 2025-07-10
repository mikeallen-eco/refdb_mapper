compile_loso_ghostblaster_outcomes <- function(df) {

# compile results for BLAST  
  df_b <- df %>%
    mutate(thresh97 = case_when((
      max_pident_nei_no_ecto >= 97 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 97 &
        top_match_nei_no_ecto == true_ncbi_name
    ) ~ "correct",
    (
      max_pident_nei_no_ecto >= 97 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 97 &
        top_match_nei_no_ecto != true_ncbi_name
    ) ~ "incorrect",!(
      max_pident_nei_no_ecto >= 97 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 97
    ) ~ "abstain"
    )) %>%
    mutate(thresh98 = case_when((
      max_pident_nei_no_ecto >= 98 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 98 &
        top_match_nei_no_ecto == true_ncbi_name) ~ "correct",
    (
      max_pident_nei_no_ecto >= 98 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 98 &
        top_match_nei_no_ecto != true_ncbi_name) ~ "incorrect",!(
      max_pident_nei_no_ecto >= 98 &
        max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto < 98) ~ "abstain"
    )) %>%
    mutate(thresh99 = case_when((
      max_pident_nei_no_ecto >= 99 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 99 &
        top_match_nei_no_ecto == true_ncbi_name) ~ "correct",
    (
      max_pident_nei_no_ecto >= 99 &
        (max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto) < 99 &
        top_match_nei_no_ecto != true_ncbi_name) ~ "incorrect",!(
      max_pident_nei_no_ecto >= 99 &
        max_pident_nei_no_ecto - max_pident_gap_nei_no_ecto < 99) ~ "abstain"
    )) %>%
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
                              score_nei_no_ecto < 90 ~ "abstain"))
  
  df_b$thresh97_i <- 1*(df_b$thresh97 %in% "incorrect")
  df_b$thresh98_i <- 1*(df_b$thresh98 %in% "incorrect")
  df_b$thresh99_i <- 1*(df_b$thresh99 %in% "incorrect")
  df_b$ecotag_i <- 1*(df_b$ecotag %in% "incorrect")
  df_b$conf90_i <- 1*(df_b$conf90 %in% "incorrect")
  df_b$thresh97_c <- 1*(df_b$thresh97 %in% "correct")
  df_b$thresh98_c <- 1*(df_b$thresh98 %in% "correct")
  df_b$thresh99_c <- 1*(df_b$thresh99 %in% "correct")
  df_b$ecotag_c <- 1*(df_b$ecotag %in% "correct")
  df_b$conf90_c <- 1*(df_b$conf90 %in% "correct")
  df_b$thresh97_a <- 1*(df_b$thresh97 %in% "abstain")
  df_b$thresh98_a <- 1*(df_b$thresh98 %in% "abstain")
  df_b$thresh99_a <- 1*(df_b$thresh99 %in% "abstain")
  df_b$ecotag_a <- 1*(df_b$ecotag %in% "abstain")
  df_b$conf90_a <- 1*(df_b$conf90 %in% "abstain")
  
  # compile results for GhostBLASTer
  df_gb <- df %>%
    mutate(thresh97 = case_when((
      max_pident_nei_ecto >= 97 &
        (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 97 &
        top_match_nei_ecto == true_ncbi_name
    ) ~ "correct",
    (
      max_pident_nei_ecto >= 97 &
        (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 97 &
        top_match_nei_ecto != true_ncbi_name
    ) ~ "incorrect",!(
      max_pident_nei_ecto >= 97 &
        (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 97
    ) ~ "abstain"
    )) %>%
    mutate(thresh98 = case_when((
      max_pident_nei_ecto >= 98 &
        (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 98 &
        top_match_nei_ecto == true_ncbi_name) ~ "correct",
      (
        max_pident_nei_ecto >= 98 &
          (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 98 &
          top_match_nei_ecto != true_ncbi_name) ~ "incorrect",!(
            max_pident_nei_ecto >= 98 &
              max_pident_nei_ecto - max_pident_gap_nei_ecto < 98) ~ "abstain"
    )) %>%
    mutate(thresh99 = case_when((
      max_pident_nei_ecto >= 99 &
        (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 99 &
        top_match_nei_ecto == true_ncbi_name) ~ "correct",
      (
        max_pident_nei_ecto >= 99 &
          (max_pident_nei_ecto - max_pident_gap_nei_ecto) < 99 &
          top_match_nei_ecto != true_ncbi_name) ~ "incorrect",!(
            max_pident_nei_ecto >= 99 &
              max_pident_nei_ecto - max_pident_gap_nei_ecto < 99) ~ "abstain"
    )) %>%
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
                              score_nei_ecto < 90 ~ "abstain"))
  
  df_gb$thresh97_i <- 1*(df_gb$thresh97 %in% "incorrect")
  df_gb$thresh98_i <- 1*(df_gb$thresh98 %in% "incorrect")
  df_gb$thresh99_i <- 1*(df_gb$thresh99 %in% "incorrect")
  df_gb$ecotag_i <- 1*(df_gb$ecotag %in% "incorrect")
  df_gb$conf90_i <- 1*(df_gb$conf90 %in% "incorrect")
  df_gb$thresh97_c <- 1*(df_gb$thresh97 %in% "correct")
  df_gb$thresh98_c <- 1*(df_gb$thresh98 %in% "correct")
  df_gb$thresh99_c <- 1*(df_gb$thresh99 %in% "correct")
  df_gb$ecotag_c <- 1*(df_gb$ecotag %in% "correct")
  df_gb$conf90_c <- 1*(df_gb$conf90 %in% "correct")
  df_gb$thresh97_a <- 1*(df_gb$thresh97 %in% "abstain")
  df_gb$thresh98_a <- 1*(df_gb$thresh98 %in% "abstain")
  df_gb$thresh99_a <- 1*(df_gb$thresh99 %in% "abstain")
  df_gb$ecotag_a <- 1*(df_gb$ecotag %in% "abstain")
  df_gb$conf90_a <- 1*(df_gb$conf90 %in% "abstain")
  
  return(list(b = df_b,
              gb = df_gb))
  
}