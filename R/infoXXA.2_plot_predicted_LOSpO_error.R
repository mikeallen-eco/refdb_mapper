plot_predicted_LOSpO_error <- function(preds_lospo){

pred_lospo_plot_i <- ggplot(preds_lospo$i$pred_df) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_x_continuous(limits = c(0,30)) +
  # scale_y_continuous(limits = c(0,50)) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "Predicted % misclassified",
       title = "Sequence from novel species\n(leave-one-species-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_lospo_plot_a <- ggplot(preds_lospo$a$pred_df) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_x_continuous(limits = c(0,30)) +
  # scale_y_continuous(limits = c(0,100)) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "Predicted % unclassified",
       title = "Sequence from novel species\n(leave-one-species-out)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

return(list(i = pred_lospo_plot_i,
            a = pred_lospo_plot_a))

}
