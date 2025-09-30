# compile metrics to choose the best combinations of markers
  # in terms of predicted p(correct) and % coverage (% species with > 0 sequences)

library(dplyr)
library(combinat)

pick_best_markers <- function(hybas_data = pf_data){
  
  # pull out performance matrices
  pc_mat <- hybas_data %>%
    select(ends_with("_p_c")) %>%
    filter(!if_any(everything(), is.na)) %>%
    as.matrix()
  
  pi_mat <- hybas_data %>%
    select(ends_with("_p_i")) %>%
    filter(!if_any(everything(), is.na)) %>%
    as.matrix()
  
  marker_names <- colnames(pc_mat)  # same set of markers for _p_c and _p_i
  
  # function to evaluate a given set of markers
  evaluate_combo_stats <- function(markers) {
    idx <- match(markers, marker_names)
    
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
      markers = sapply(combos, paste, collapse = "+")
    ) %>%
      bind_cols(bind_rows(stats_list)) %>%
      arrange(desc(median_p_correct))
  }
  
  best_df <- best_combos(1) %>% 
    bind_rows(best_combos(2),
              best_combos(3)) %>%
    arrange(desc(median_p_correct)) %>%
    mutate(markers = gsub("_p_c", "", markers),
           markers = gsub("\\+", " + ", markers))
  
  return(best_df)
}

# library(dplyr)
# library(combinat)   # or use utils::combn
# 
# pick_best_markers <- function(hybas_data = pf_data){
#   
#   pc_mat <- hybas_data %>%
#     select(ends_with("_p_c")) %>%
#     filter(!if_any(everything(), is.na)) %>%
#     as.matrix()
#   
#   marker_names <- colnames(pc_mat)
#   
#   # function to evaluate a given set of markers
#   evaluate_combo_stats <- function(markers) {
#     idx <- match(markers, marker_names)
#     combo_scores <- apply(pc_mat[, idx, drop = FALSE], 1, max, na.rm = TRUE)
#     tibble(
#       median_p_correct = median(combo_scores, na.rm = TRUE),
#       q25_p_correct    = quantile(combo_scores, 0.25, na.rm = TRUE),
#       q75_p_correct    = quantile(combo_scores, 0.75, na.rm = TRUE),
#       pct_cov          = mean(combo_scores > 0, na.rm = TRUE)
#     )
#   }
#   
#   best_combos <- function(pc_mat, N) {
#     combos <- combn(marker_names, N, simplify = FALSE)
#     
#     stats_list <- lapply(combos, evaluate_combo_stats)
#     
#     tibble(
#       markers = sapply(combos, paste, collapse = "+")
#     ) %>%
#       bind_cols(bind_rows(stats_list)) %>%
#       arrange(desc(median_p_correct))
#   }
#   
#   best_df <- best_combos(pc_mat, 1) %>% 
#     bind_rows(best_combos(pc_mat, 2),
#               best_combos(pc_mat, 3)) %>%
#     arrange(desc(median_p_correct)) %>%
#     mutate(markers = gsub("_p_c", "", markers),
#            markers = gsub("\\+", " + ", markers))
#   
#   return(best_df)
# }


# pick_best_markers <- function(hybas_data = pf_data){
#   
#   pc_mat <- hybas_data %>%
#     select(ends_with("_p_c")) %>%
#     filter(!if_any(everything(), is.na)) %>%
#     as.matrix()
#   
#   # try all marker combinations
#   
#   marker_names <- colnames(pc_mat)
#   
#   # function to evaluate a given set of markers
#   evaluate_combo <- function(markers) {
#     idx <- match(markers, marker_names)
#     combo_scores <- apply(pc_mat[, idx, drop = FALSE], 1, max, na.rm = TRUE)
#     median(combo_scores, na.rm = TRUE)  # overall average
#   }
#   
#   evaluate_combo_pct <- function(markers) {
#     idx <- match(markers, marker_names)
#     combo_scores <- apply(pc_mat[, idx, drop = FALSE], 1, max, na.rm = TRUE)
#     mean(combo_scores > 0, na.rm = TRUE)
#   }
#   
#   best_combos <- function(pc_mat, N) {
#     marker_names <- colnames(pc_mat)
#     combos <- combn(marker_names, N, simplify = FALSE)
#     scores <- sapply(combos, evaluate_combo)
#     p_c_df <- tibble(
#       markers = sapply(combos, paste, collapse = "+"),
#       median_p_correct = scores
#     ) %>%
#       arrange(desc(median_p_correct))
#     
#     scores_pct <- sapply(combos, evaluate_combo_pct)
#     pct_df <- tibble(
#       markers = sapply(combos, paste, collapse = "+"),
#       pct_cov = scores_pct
#     ) %>%
#       arrange(desc(pct_cov))
#     
#     final_df <- p_c_df %>%
#       left_join(pct_df, 
#                 by = join_by(markers))
#     
#   }
#   
#   best_df <- best_combos(pc_mat, 1) %>% 
#     bind_rows(best_combos(pc_mat, 2),
#               best_combos(pc_mat, 3)) %>%
#     arrange(desc(median_p_correct)) %>%
#     mutate(markers = gsub("_p_c", "", markers)) %>%
#     mutate(markers = gsub("\\+", " + ", markers))
#   
#   return(best_df)
#   
# }