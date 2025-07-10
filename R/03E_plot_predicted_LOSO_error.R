plot_predicted_LOSO_error <- function(preds){

library(ggplot2)
pred_loso_plot <- ggplot(preds$pred_df) +
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

# ggsave("figures/LOSO_predicted_%_misclassified_BLAST99.png",
#        height = 4, width = 6, dpi = 400, bg = "white")

return(pred_loso_plot)

}