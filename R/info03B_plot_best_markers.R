# plot performance of marker combinations using the df generated in pick_best_markers()

library(ggplot2)
library(forcats)

plot_best_markers <- function(best_markers_df = pf_best, 
                              metric = "c"){
  
  if(metric == "c"){
    
  best_rubrics <- best_markers_df %>%
    group_by()
  best_markers_df %>%
    mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
    arrange(desc(median_p_correct), num_markers) %>%
    mutate(order = length(num_markers):1) %>%
    mutate(markers = forcats::fct_reorder(markers, order)) %>%
    ggplot() +
    geom_errorbar(aes(y = markers, xmin = 100*q25_p_correct, 
                      xmax = 100*q75_p_correct), width = 0, linewidth = 1) +
    geom_point(aes(y = markers, x = 100*median_p_correct), size = 6) +
    geom_point(aes(y = markers, x = 100*pct_cov), size = 6, color = "steelblue") +
    theme_bw() +
    scale_x_continuous(limits = c(0,100), breaks = seq(0,100, by = 10)) +
    theme(text = element_text(size = 16)) +
    labs(x = "Median p(assigned ∩ correct) or \n total species coverage (%)",
         y = "",
         title = "Marker combinations")
  }
  
  if(metric == "i"){
  x_upper <- max(max(best_markers_df$q75_p_incorrect), 20)
  best_markers_df %>%
    mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
    arrange(desc(median_p_incorrect), num_markers) %>%
    mutate(order = 1:length(num_markers)) %>%
    mutate(markers = forcats::fct_reorder(markers, order)) %>%
    ggplot() +
    geom_errorbar(aes(y = markers, xmin = 100*q25_p_incorrect, 
                      xmax = 100*q75_p_incorrect), width = 0, linewidth = 1) +
    geom_point(aes(y = markers, x = 100*median_p_incorrect), size = 6) +
    theme_bw() +
    scale_x_continuous(limits = c(0,x_upper), breaks = seq(0,100, by = 10)) +
    theme(text = element_text(size = 16)) +
    labs(x = "Median p(assigned ∩ incorrect)",
         y = "",
         title = "Marker combinations")    
  }
  
}
