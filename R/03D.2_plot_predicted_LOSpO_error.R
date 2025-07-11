plot_predicted_LOSpO_error <- function(preds){

library(ggplot2)
pred_lospo_plot_i <- ggplot(preds$i$pred_df) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_x_continuous(limits = c(0,30)) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "Predicted % misclassified",
       title = "Leave-one-species-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

pred_lospo_plot_i <- ggplot(preds$i$pred_df) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_x_continuous(limits = c(0,30)) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "Predicted % unclassified",
       title = "Leave-one-species-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

return(list(i = pred_lospo_plot_i,
            a = pred_lospo_plot_a))

}