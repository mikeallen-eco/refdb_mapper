plot_predicted_LOSO_error <- function(preds_loso){

pred_loso_plot_i <- preds_loso$i$pred_df %>%
  mutate(preds = case_when(preds > 0.1 ~ 0.1,
                           TRUE ~ preds)) %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno", limits = c(0,11), 
                       breaks = c(0,2,4,6,8,10), labels = c(0,2,4,6,8,"10+")) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(1,10), labels = 1:10, breaks = 1:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified",
       title = "Novel sequence\n(leave-one-sequence-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_loso_plot_a <- preds_loso$a$pred_df %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(1,10), labels = 1:10, breaks = 1:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nunclassified",
       title = "Novel sequence\n(leave-one-sequence-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_loso_plot_c <- preds_loso$c$pred_df %>%
  ggplot() +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(1,10), labels = 1:10, breaks = 1:10) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\ncorrect",
       title = "Novel sequence\n(leave-one-sequence-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

return(list(i = pred_loso_plot_i,
            a = pred_loso_plot_a,
            c = pred_loso_plot_c))

}
