# ---- Step 3: plot % correct vs. % incorrect

library(ggplot2)
library(wesanderson)

scatterplot_best_markers <- function(best_df = ct_best,
                                     top_n = 1,
                                     num_markers = c(1,2,3)){
  
  if(is.null(top_n)){top_n <- 2000}
  
  num_markers_vector <- num_markers
  
  plot_data_best <- best_df %>%
    filter(num_markers %in% num_markers_vector) %>%
    group_by(as.factor(num_markers)) %>%
    arrange(desc(std_score)) %>%
    slice_head(n = top_n) %>% 
    distinct()
  
  plot_data_best_labels <- best_df %>%
    filter(num_markers %in% num_markers_vector) %>%
    group_by(as.factor(num_markers)) %>%
    arrange(desc(std_score)) %>%
    slice_head(n = 1) %>% 
    distinct()
  
  # Compute global limits once
  score_limits <- range(best_df$std_score, na.rm = TRUE)
  
  # Extract top 3 marker texts
  marker_texts <- plot_data_best_labels$markers[1:3]
  rubric_texts <- plot_data_best_labels$rubric[1:3]
  
  # Build the multiline annotation
  annotation_text <- paste0(
    "One marker: ", marker_texts[1], "( ", rubric_texts[1], ")\n",
    "Two markers: ", marker_texts[2], " (", rubric_texts[2], ")\n",
    "Three markers: ", marker_texts[3], " (", rubric_texts[2], ")"
  )
  
  # Add to your existing plot code
  plot <- plot_data_best %>%
    ggplot() +
    geom_point(aes(x = median_p_correct, y = median_p_incorrect, 
                   color = std_score, size = num_markers), 
               data = best_df, alpha = 0.65) +
    geom_errorbar(aes(ymin = q25_p_incorrect, ymax = q75_p_incorrect,
                      x = median_p_correct, y = median_p_incorrect),
                  width = 0, linewidth = 0.75, alpha = 0.65, color = "black") +
    geom_errorbar(aes(xmin = q25_p_correct, xmax = q75_p_correct,
                      x = median_p_correct, y = median_p_incorrect),
                  width = 0, linewidth = 0.75, alpha = 0.65,
                  orientation = "y", color = "black") +
    geom_point(aes(x = median_p_correct, y = median_p_incorrect, 
                   size = num_markers, fill = std_score),
               shape = 21, color = "black") +
    scale_color_viridis_c(name = "Std. score", limits = score_limits) +
    scale_fill_viridis_c(name = "Std. score", limits = score_limits, guide = "none") +
    scale_size_continuous(
      breaks = 1:3,
      limits = c(1, 3),
      range = c(3, 8)
    ) +
    labs(
      x = "P (assigned ∩ correct)",
      y = "P (assigned ∩ incorrect)",
      size = "No.\nmarkers"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14)) +
    # Add annotation text to top right corner
    annotate("text",
             x = Inf, y = Inf,          # top-right corner
             label = annotation_text,
             hjust = 1.1, vjust = 1.1,  # adjust padding inward
             size = 3.5, lineheight = 1.1,
             fontface = "italic")
  
  suppressWarnings(plot)
  
}

