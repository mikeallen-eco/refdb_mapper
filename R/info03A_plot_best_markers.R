# plot performance of marker combinations using the df generated in pick_best_markers()

library(ggplot2)
library(forcats)

plot_best_markers <- function(best_markers_df = pf_best, 
                              metric = "c",
                              rubrics = c("blast97", "blast98", 
                                          "blast99", "ecotag",
                                          "rdp70", "rdp80",
                                          "rdp90", "rdp95")){
  
  if(metric == "c"){
    
    bestrubrics <- best_markers_df %>%
      filter(rubric %in% rubrics) %>%
      mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
      arrange(rubric, desc(median_p_correct), num_markers) %>%
      group_by(rubric) %>%
      slice_head(n = 1) %>%
      arrange(desc(median_p_correct)) %>%
      ungroup() %>%
      slice_head(n = 2)
    
    top_rubric <- bestrubrics$rubric[1]
    
    # Ordering of markers just for the top rubric
    marker_order <- best_markers_df %>%
      filter(rubric == top_rubric) %>%
      mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
      arrange(desc(median_p_correct), num_markers) %>%
      pull(markers)
    
    plot_df <- best_markers_df %>%
      filter(rubric %in% bestrubrics$rubric) %>%
      mutate(rubric = factor(rubric, levels = bestrubrics$rubric)) %>%
      mutate(markers = factor(markers, levels = rev(marker_order)))   # top-to-bottom
    
    plot <- plot_df %>%
      ggplot() +
      geom_errorbar(
        aes(y = markers,
            xmin = 100 * q25_p_correct,
            xmax = 100 * q75_p_correct,
            group = rubric, color = rubric),
        position = position_dodge(width = 1),
        width = 0, linewidth = 1
      ) +
      geom_point(
        aes(y = markers, x = 100 * median_p_correct, 
            shape = rubric, color = rubric),
        position = position_dodge(width = 1),
        size = 6
      ) +
      geom_point(
        aes(y = markers, x = 100 * pct_cov),
        size = 6, color = "steelblue", shape = "+"
      ) +
      theme_bw() +
      scale_color_manual(values = c("firebrick", "black")) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
      theme(text = element_text(size = 16)) +
      labs(
        x = "Median p(assigned ∩ correct) or \n + = species coverage (%)",
        y = "",
        title = "Marker combinations",
        shape = ""
      ) +
      guides(color = "none")
    
  }
  
  if(metric == "i"){
    
    bestrubrics <- best_markers_df %>%
      mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
      group_by(rubric) %>%
      # get each rubric’s best marker (lowest incorrect)
      arrange(median_p_incorrect, num_markers, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(median_p_incorrect) %>%
      slice_head(n = 2)   # keep two best rubrics
    
    lowest_rubric <- bestrubrics$rubric[1]
    
    marker_order <- best_markers_df %>%
      filter(rubric == lowest_rubric) %>%
      mutate(num_markers = stringr::str_count(markers, "\\+") + 1) %>%
      arrange(median_p_incorrect, num_markers) %>%
      pull(markers)
    
    x_upper <- max(max(bestrubrics$q75_p_incorrect), 20)
    
    plot_df <- best_markers_df %>%
      filter(rubric %in% bestrubrics$rubric) %>%
      mutate(rubric = factor(rubric, levels = bestrubrics$rubric)) %>%
      mutate(markers = factor(markers, levels = rev(marker_order)))  # top = lowest incorrect
    
    plot <- plot_df %>%
      ggplot() +
      geom_errorbar(
        aes(y = markers,
            xmin = 100 * q25_p_incorrect,
            xmax = 100 * q75_p_incorrect,
            group = rubric, color = rubric),
        position = position_dodge(width = 1),
        width = 0, linewidth = 1
      ) +
      geom_point(
        aes(y = markers, x = 100 * median_p_incorrect, 
            shape = rubric, color = rubric),
        position = position_dodge(width = 1),
        size = 6
      ) +
      theme_bw() +
      scale_color_manual(values = c("firebrick", "black")) +
      scale_x_continuous(limits = c(0, x_upper), breaks = seq(0, 100, by = 10)) +
      theme(text = element_text(size = 16)) +
      labs(
        x = "Median p(assigned ∩ incorrect)",
        y = "",
        title = "Marker combinations",
        shape = ""
      ) +
      guides(color = "none")
 
  }
  
  return(plot)
}
