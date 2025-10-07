plot_predicted_error_both <- function(preds_loso, 
                                      preds_lospo){

  
## --- format p(incorrect) data for plotting
  
  pred_loso_plot_data_i <- preds_loso$i$pred_df
      
  pred_lospo_plot_data_i <- expand.grid(nnd = unique(preds_lospo$i$pred_df$nnd),
                           n_seqs = seq(0, 1, by = 0.090909)[2:11]) %>%
        mutate(preds = rep(unique(preds_lospo$i$pred_df$preds), 10))
  
  pred_both_plot_data_i <- pred_loso_plot_data_i %>%
    bind_rows(pred_lospo_plot_data_i) %>%
    arrange(n_seqs, nnd, preds) %>%
    mutate(preds = case_when(preds > 0.1 ~ 0.1, TRUE ~ preds))
  
## --- format p(abstain) data for plotting
  
  pred_loso_plot_data_a <- preds_loso$a$pred_df
  
  pred_lospo_plot_data_a <- expand.grid(nnd = unique(preds_lospo$a$pred_df$nnd),
                                        n_seqs = seq(0, 1, by = 0.090909)[2:11]) %>%
    mutate(preds = rep(unique(preds_lospo$a$pred_df$preds), 10))
  
  pred_both_plot_data_a <- pred_loso_plot_data_a %>%
    bind_rows(pred_lospo_plot_data_a) %>%
    arrange(n_seqs, nnd, preds)

## --- format p(correct) data for plotting
  
  pred_loso_plot_data_c <- preds_loso$c$pred_df
  
  pred_lospo_plot_data_c <- expand.grid(nnd = unique(preds_loso$c$pred_df$nnd),
                                        n_seqs = seq(0, 1, by = 0.090909)[2:11]) %>%
    mutate(preds = rep(0, 1000))
  
  pred_both_plot_data_c <- pred_loso_plot_data_c %>%
    bind_rows(pred_lospo_plot_data_c) %>%
    arrange(n_seqs, nnd, preds) %>%
    mutate(preds = case_when(preds < 0.5 ~ 0.5, TRUE ~ preds))
  
## --- Plot p(incorrect) including both loso and lospo models
(pred_both_plot_i <- pred_both_plot_data_i %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno", limits = c(0,10), 
                       breaks = c(0,2,4,6,8,10), labels = c(0,2,4,6,8,"10+")) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,10), labels = 0:10, breaks = 0:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "%\nassigned &\nincorrect",
       title = "Novel sequence or species\n(leave-one-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title.align = 0.5))

(pred_both_plot_a <- pred_both_plot_data_a %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,10), labels = 0:10, breaks = 0:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "%\nunassigned",
       title = "Novel sequence or species\n(leave-one-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title.align = 0.5))

(pred_both_plot_c <- pred_both_plot_data_c %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  # scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  scale_fill_viridis_c(option = "inferno", limits = c(50,100), 
                         breaks = c(50,60,70,80,90,100), labels = c("≤50",60,70,80,90,100)) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,10), labels = 0:10, breaks = 0:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "%\nassigned &\ncorrect",
       title = "Novel sequence or species\n(leave-one-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title.align = 0.5))

return(list(i = pred_both_plot_i,
            a = pred_both_plot_a,
            c = pred_both_plot_c))

}
