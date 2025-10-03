# ---- Step 3: plot % correct vs. % incorrect

library(ggplot2)
library(wesanderson)

scatterplot_best_markers <- function(best_df = ct_best,
                                     top_n = NA,
                                     num_markers = c(1,2)){
  
  if(is.null(top_n)){top_n <- 2000}
  
  num_markers_vector <- num_markers
  
  plot_data_best <- best_df %>%
    mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
    filter(num_markers %in% num_markers_vector)
  
  plot_data_best_c <- plot_data_best %>%
    arrange(desc(median_p_correct)) %>%
    slice_head(n = top_n)
  
  plot_data_best_i <- plot_data_best %>%
    arrange(median_p_incorrect) %>%
    slice_head(n = top_n)
  
  plot_data <- plot_data_best_c %>%
    bind_rows(plot_data_best_i) %>%
    distinct()
  
  plot <- plot_data %>%
    ggplot() +
    geom_errorbar(aes(ymin = q25_p_incorrect, ymax = q75_p_incorrect,
                      x = median_p_correct, y = median_p_incorrect, color = rubric),
                  width = 0, linewidth = 0.75, alpha = 0.65) +
    geom_errorbarh(aes(xmin = q25_p_correct, xmax = q75_p_correct,
                       x = median_p_correct, y = median_p_incorrect, color = rubric),
                   width = 0, linewidth = 0.75, alpha = 0.65) +
    geom_point(aes(x = median_p_correct, y = median_p_incorrect, 
                   color = rubric, size = as.factor(num_markers)), 
               shape = 21, fill = "white") +
    scale_color_manual(values = c(wes_palettes$Darjeeling1, wes_palettes$Darjeeling2)) +
    labs(x = "P (assigned ∩ correct)",
         y = "P (assigned ∩ incorrect)",
         color = "Rubric", size = "No.\nmarkers") +
    theme_minimal() +
    theme(text = element_text(size = 14))
  
  if(max_num_markers > 1){ 
    plot <- plot  +
      geom_errorbar(aes(ymin = q25_p_incorrect, ymax = q75_p_incorrect,
                        x = median_p_correct, y = median_p_incorrect, color = rubric),
                    width = 0, linewidth = 0.75, alpha = 0.65, 
                    data = plot_data %>% filter(num_markers %in% 2)) +
      geom_errorbarh(aes(xmin = q25_p_correct, xmax = q75_p_correct,
                         x = median_p_correct, y = median_p_incorrect, color = rubric),
                     width = 0, linewidth = 0.75, alpha = 0.65, 
                     data = plot_data %>% filter(num_markers %in% 2)) +
      geom_point(aes(x = median_p_correct, y = median_p_incorrect, 
                     color = rubric, size = as.factor(num_markers)), 
                 shape = 21, fill = "white", data = plot_data %>% filter(num_markers %in% 2)) +
      geom_errorbar(aes(ymin = q25_p_incorrect, ymax = q75_p_incorrect,
                        x = median_p_correct, y = median_p_incorrect, color = rubric),
                    width = 0, linewidth = 0.75, alpha = 0.65, 
                    data = plot_data %>% filter(num_markers %in% 1)) +
      geom_errorbarh(aes(xmin = q25_p_correct, xmax = q75_p_correct,
                         x = median_p_correct, y = median_p_incorrect, color = rubric),
                     width = 0, linewidth = 0.75, alpha = 0.65, 
                     data = plot_data %>% filter(num_markers %in% 1)) +
      geom_point(aes(x = median_p_correct, y = median_p_incorrect, 
                     color = rubric, size = as.factor(num_markers)), 
                 shape = 21, fill = "white", data = plot_data %>% filter(num_markers %in% 1))
  }
  
  suppressWarnings(plot)
  
}
