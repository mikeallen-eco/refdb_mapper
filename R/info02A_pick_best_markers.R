# compile metrics to choose the best combinations of markers for each rubric
  # in terms of predicted p(correct) and % coverage (% species with > 0 sequences)

library(dplyr)
library(combinat)

pick_best_markers <- function(hybas_data = nj_data,
                              rubrics = c("blast97", "blast98", "blast99", "ecotag",
                                          "rdp70", "rdp80", "rdp90", "rdp95")){
  
  rubric_list <- lapply(seq_along(rubrics), function(i){
    
    valid_rows <- hybas_data %>%
      select(contains(paste0(rubrics[i], "_p_c")),
             contains(paste0(rubrics[i], "_p_i"))) %>%
      drop_na()
    
    pc_mat <- valid_rows %>%
      select(contains(paste0(rubrics[i], "_p_c"))) %>%
      as.matrix()
    
    pi_mat <- valid_rows %>%
      select(contains(paste0(rubrics[i], "_p_i"))) %>%
      as.matrix()
    
    
  # # pull out performance matrices
  # pc_mat <- hybas_data %>%
  #   select(contains(paste0(rubrics[i], "_p_c"))) %>%
  #   filter(!if_any(everything(), is.na)) %>%
  #   as.matrix()
  # 
  # pi_mat <- hybas_data %>%
  #   select(contains(paste0(rubrics[i], "_p_i"))) %>%
  #   filter(!if_any(everything(), is.na)) %>%
  #   as.matrix()
  
  marker_names <- colnames(pc_mat)  # same set of markers for _p_c and _p_i
  
  # function to evaluate a given set of markers
  evaluate_combo_stats <- function(marker_combos) {
    idx <- match(marker_combos, marker_names)
    
    # combine by taking max across markers for p_c and p_i
    combo_scores_c <- apply(pc_mat[, idx, drop = FALSE], 1, max, na.rm = TRUE)
    combo_scores_i <- apply(pi_mat[, idx, drop = FALSE], 1, min, na.rm = TRUE)
    
    tibble(
      median_p_correct = median(combo_scores_c, na.rm = TRUE),
      q25_p_correct    = quantile(combo_scores_c, 0.25, na.rm = TRUE),
      q75_p_correct    = quantile(combo_scores_c, 0.75, na.rm = TRUE),
      pct_cov          = mean(combo_scores_c > 0, na.rm = TRUE),
      median_p_incorrect = median(combo_scores_i, na.rm = TRUE),
      q25_p_incorrect    = quantile(combo_scores_i, 0.25, na.rm = TRUE),
      q75_p_incorrect    = quantile(combo_scores_i, 0.75, na.rm = TRUE)
    )
  }
  
  best_combos <- function(N) {
    combos <- combn(marker_names, N, simplify = FALSE)
    stats_list <- lapply(combos, evaluate_combo_stats)
    
    tibble(
      markers = sapply(combos, paste, collapse = " + ")
    ) %>%
      bind_cols(bind_rows(stats_list)) %>%
      arrange(desc(median_p_correct))
  }
  
  best_df <- best_combos(1) %>% 
    bind_rows(best_combos(2),
              best_combos(3)) %>%
    arrange(desc(median_p_correct)) %>%
    mutate(markers = gsub(paste0("_", rubrics[i], "_p_c"), "", markers),
           rubric = rubrics[i])
  
  return(best_df)
  })
  
  names(rubric_list) <- rubrics
  
  final_df <- bind_rows(rubric_list)
  
  return(final_df)
}
