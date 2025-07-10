plot_predicted_LOSpO_error <- function(preds){

library(ggplot2)
pred_lospo_plot <- ggplot(preds$pred_df) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,25)) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "Predicted % misclassified",
       title = "Leave-one-species-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

return(pred_lospo_plot)

}