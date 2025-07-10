plot_predicted_error <- function(preds){

library(ggplot2)
thresh99_i <- ggplot(preds$pred_df_thresh99_i) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:10, breaks = 2:11) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified\n(BLAST 99%\nthreshold)",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# ggsave("figures/LOSO_predicted_%_misclassified_BLAST99.png",
#        height = 4, width = 6, dpi = 400, bg = "white")

ecotag_i <- ggplot(preds$pred_df_ecotag_i) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:11, breaks = 2:12) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified\n(ecotag)",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# ggsave("figures/LOSO_predicted_%_misclassified_ecotag.png", 
#        height = 4, width = 6, dpi = 400, bg = "white")

conf90_blast <- ggplot(preds$pred_df_conf90_i) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:11, breaks = 2:12) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified\n(BLAST\n90%)",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# ggsave("figures/LOSO_predicted_%_misclassified_ghostblaster90.png", 
#        height = 4, width = 6, dpi = 400, bg = "white")

conf90_ghostblaster <- ggplot(preds$pred_df_conf90_gb_i) +
  geom_tile(aes(x = nnd, y = n_seqs, fill = 100*preds)) +
  scale_fill_viridis_c(option = "inferno") +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(2,11), labels = 1:11, breaks = 2:12) +
  labs(x = "Nearest evolutionary neighbor (MY)",
       y = "No. sequences in reference database",
       fill = "Predicted %\nmisclassified\n(GhostBLASTer\n90%)",
       title = "Leave-one-sequence-out") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# ggsave("figures/LOSO_predicted_%_misclassified_ghostblaster90.png", 
#        height = 4, width = 6, dpi = 400, bg = "white")

return(list(thresh99_i = thresh99_i,
            ecotag_i = ecotag_i,
            conf90_blast = conf90_blast,
            conf90_ghostblaster = conf90_ghostblaster))

}