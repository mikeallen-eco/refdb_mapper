plot_predicted_LOSO_error <- function(preds){

library(ggplot2)
pred_loso_plot_i <- ggplot(preds$i$pred_df) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:10, breaks = 2:11) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_loso_plot_a <- ggplot(preds$a$pred_df) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:10, breaks = 2:11) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nunclassified",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_loso_plot_c <- ggplot(preds$c$pred_df) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:10, breaks = 2:11) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\ncorrect",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

return(list(pred_loso_plot_i,
            pred_loso_plot_a,
            pred_loso_plot_c))

}