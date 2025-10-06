# plot performance of marker combinations using the df generated in pick_best_markers()

library(ggplot2)
library(forcats)

plot_best_markers <- function(best_markers_df = ct_best,
                              num_markers = c(1,2,3),
                              top_n_rubrics = 2,
                              rubrics = c("blast97", "blast98", 
                                          "blast99", "ecotag",
                                          "rdp70", "rdp80",
                                          "rdp90", "rdp95")){
  
  num_markers_vector <- num_markers
  
    bestrubrics <- best_markers_df %>%
      filter(rubric %in% rubrics) %>%
      filter(num_markers %in% num_markers_vector) %>%
      arrange(rubric, desc(std_score), num_markers) %>%
      group_by(rubric) %>%
      slice_head(n = 1) %>%
      arrange(desc(std_score)) %>%
      ungroup() %>%
      slice_head(n = top_n_rubrics)
    
    top_rubric <- bestrubrics$rubric[1]
    
    # Ordering of markers just for the top rubric
    marker_order <- best_markers_df %>%
      filter(rubric == top_rubric) %>%
      filter(num_markers %in% num_markers_vector) %>%
      arrange(desc(std_score), num_markers) %>%
      pull(markers)
    
    plot_df_cor <- best_markers_df %>%
      filter(rubric %in% bestrubrics$rubric) %>%
      filter(num_markers %in% num_markers_vector) %>%
      select(markers, rubric, median = median_p_correct, q25 = q25_p_correct, 
             q75 = q75_p_correct, pct_cov) %>%
      mutate(metric = "p(assigned ∩ correct)")
    
    plot_df_incor <- best_markers_df %>%
      filter(rubric %in% bestrubrics$rubric) %>%
      filter(num_markers %in% num_markers_vector) %>%
      select(markers, rubric, median = median_p_incorrect, q25 = q25_p_incorrect, 
             q75 = q75_p_incorrect) %>%
      mutate(pct_cov = NA,
             metric = "p(assigned ∩ incorrect)")
      
    plot_df <- plot_df_cor %>%
      bind_rows(plot_df_incor) %>%
      mutate(rubric = factor(rubric, levels = bestrubrics$rubric)) %>%
      mutate(markers = factor(markers, levels = rev(marker_order)))   # top-to-bottom
    
    plot <- plot_df %>%
      ggplot() +
      geom_errorbar(
        aes(y = markers,
            xmin = 100 * q25,
            xmax = 100 * q75,
            group = rubric, color = rubric),
        position = position_dodge(width = 1),
        width = 0, linewidth = 1
      ) +
      geom_point(
        aes(y = markers, x = 100 * median, 
            shape = rubric, color = rubric),
        position = position_dodge(width = 1),
        size = 6
      ) +
      geom_point(
        aes(y = markers, x = 100 * pct_cov),
        size = 6, color = "steelblue", shape = "+"
      ) +
      facet_wrap(~metric, scales = "free_x") +
      theme_bw() +
      # scale_color_manual(values = c("darkgray", "orangered")) +
      # scale_x_continuous(breaks = seq(0, 100, by = 10)) + # limits = c(0, 100), 
      theme(text = element_text(size = 16)) +
      labs(
        x = "Median probability (points) or\n% species coverage (blue +)",
        y = "",
        title = "Marker combinations",
        shape = ""
      ) +
      scale_color_manual(values = c("gray30", "orangered", "olivedrab"), name = "Rubric") +
      scale_shape_manual(values = c(16, 17, 18), name = "Rubric") #+
      # guides(
        # color = guide_legend(override.aes = list(shape = c(16, 17, 18)))
      # )
    
  return(plot)
}
