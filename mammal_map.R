# ---- SET-UP
source("R/setup.R")
source("R/settings_Vences_16S.R")

# ---- Step 0 Make taxon-specific reference database for LOSO analysis (takes ~ 10 min)

curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

# ---- Step 1 LOSO analysis (takes hours, saves csv files for convenience)

LOSO_ghostblaster(LOSO_refdb = LOO_refdb_path,
                  out = out_path, marker = marker, start_seq = 1)

LOSpO_ghostblaster(LOSpO_refdb = LOO_refdb_path,
                  out = out_path, marker = marker, start_seq = 1)

# ---- Step 2 identify & map ghosts

hydrobasin_ghosts <- identify_ghosts(LOO_refdb_path = LOO_refdb_path)

(hydrobasin_ghosts_tally <- tally_ghosts(hydrobasin_ghosts))

seq_info_summarized_by_hydrobasin <- summarize_by_hydrobasin(hydrobasin_ghosts)

# ---- Step 3 fit model of error rates based on LOSO analysis, NND, and n_seqs (takes a few min)

NND_per_sp_within_refdb <- get_NND_per_sp_within_refdb(
  refdb = LOO_refdb_path)

loo_outcomes <- get_loo_outcomes(loso_gb_path = paste0(out_path, "loso"), # get tax assign outcomes
                                  lopso_gb_path = paste0(out_path, "lospo"),
                                  refdb_nnd = NND_per_sp_within_refdb) # 2 min

preds_loso_lospo <- get_preds_loso_lospo(assign = "BLAST", # fit models
                                         accuracy_metric = "thresh98",
                                         outcomes = loo_outcomes)

# ---- Step 4 get NND, n_seqs, & error rate predictions for each sp within hydrobasins
  
# hybas_nnd <- get_NDD_per_sp_all_hydrobasins() # ~ .74 s per hydrobasin
# saveRDS(hybas_nnd, paste0(out_path, "hybas_nnd_world_fixed.rds"))
hybas_nnd <- readRDS(paste0(out_path, "hybas_nnd_world_fixed.rds"))

hybas_error_data <- format_hybas_error_data(hybas_nnd_df = hybas_nnd,
                                    ref_path = LOO_refdb_path,
                                    preds_list = preds_loso_lospo)

hybas_pred_map_sf <- make_hybas_pred_map_sf(hydrobasin_map = hydrobasin_map,
                                            hybas_pred = hybas_error_data,
                                            ghost_info = seq_info_summarized_by_hydrobasin)

# ---- Step 5 example in one hydrobasin (Delaware River)

dewa <- -1529894232; rar <- -1529894602

# geotax <- fread("data/geotax.csv")
dewa_ghosts <- hydrobasin_ghosts %>%
  filter(HYBAS_ID %in% dewa)

(dewa_ghosts_tally <- tally_ghosts(dewa_ghosts))
dewa_info_summarized <- seq_info_summarized_by_hydrobasin %>%
  filter(HYBAS_ID %in% dewa)

dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

hybas_NE <- subset_hydrobasins_states(hydrobasin_map = hybas_pred_map_sf)
# Get US states from Natural Earth
us_states <- ne_states(country = "United States of America", returnclass = "sf")

# List of Northeastern states (per US Census Bureau definition)
northeast_states <- c(
  "Connecticut", "Maine", "Massachusetts", "New Hampshire", "Virginia", "West Virginia",
  "Rhode Island", "Vermont", "New Jersey", "New York", "Pennsylvania", "Maryland", "Delaware"
)

# Filter for northeastern states
ne_us_sf <- us_states %>%
  filter(name %in% northeast_states)

ggplot(ne_us_sf) +
  geom_sf() +
  geom_sf(data = dewa_hybas_pred_map_sf, fill = "red", alpha = 0.7) +
  theme_minimal() +
  theme(text = element_text(size = 18))

ggsave("figures/conceptual/NE_map.png", height = 6, width = 6, dpi = 400)

dewa_hybas_error_data <- hybas_error_data %>%
  filter(HYBAS_ID %in% dewa)

fig_df <- dewa_hybas_error_data %>%
  select(order, geo_name, n_seqs, nnd, 
         p_misclass = preds_i, p_unclass = preds_a, p_correct = preds_c) %>%
  mutate(nnd = round(nnd),
         p_misclass = round(p_misclass, 1),
         p_unclass = round(p_unclass, 1),
         p_correct = round(p_correct, 1))

library(knitr)
library(kableExtra)

# Optional: Rename columns to be more readable
fig_df_pretty <- fig_df %>%
  arrange(desc(p_misclass)) %>%
  dplyr::rename(
    `Order` = order,
    `Species` = geo_name,
    `# Sequences` = n_seqs,
    `NND` = nnd,
    `P(Misclass)` = p_misclass,
    `P(Unclass)` = p_unclass,
    `P(Correct)` = p_correct
  )

fig_df_pretty_least <- fig_df_pretty %>%
  arrange(`P(Misclass)`)

fig_df_pretty_alpha <- fig_df_pretty %>%
  arrange(`Species`)

# Print nice table
# most error prone
most_table_html <- knitr::kable(fig_df_pretty, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# least error prone
least_table_html <- knitr::kable(fig_df_pretty_least, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# alphabetical
alpha_table_html <- knitr::kable(fig_df_pretty_alpha, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# Save to file
kableExtra::save_kable(most_table_html, "figures/conceptual/eg_most_error_table.html")
kableExtra::save_kable(least_table_html, "figures/conceptual/eg_least_error_table.html")
kableExtra::save_kable(alpha_table_html, "figures/conceptual/eg_alpha_table.html")

# ---- Step 99 make maps and plots

(ghost_plots <- map_ghosts(df = seq_info_summarized_by_hydrobasin,
                           hydrobasin_map = hydrobasin_map))

(med_num_seqs_plot <- map_ghosts(df = seq_info_summarized_by_hydrobasin,
                                 hydrobasin_map = hydrobasin_map, 
                                 to_plot = "med_seqs"))

predicted_loso_lospo_error_plots <- plot_predicted_loso_lopso_error(preds = preds_loso_lospo)

hybas_misclass_rate_maps <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf)

hybas_mean_nnd_map <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf,
                                                    to_plot = "nnd")

err_scatter <- error_rate_scatterplots(hybas_pred_map_sf)

# plot forecasted improvement if all spp had ≥ 10 seqs
hybas_misclass_rate_maps10 <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf%>%
                                                              select(-starts_with("mean_")) %>%
                                                              rename_with(~ gsub("mean10", "mean", .),
                                                                          starts_with("mean10")))

# ---- Step 100 save maps and plots to png files

save_num_mammals(ghost_plots$num_all)
save_mean_nnd(hybas_mean_nnd_map$nnd)
save_num_pct_molecular_ghosts(ghost_plots$num_ghosts/ghost_plots$pct_ghosts)
save_median_num_seqs_non_ghosts(med_seqs$med_seqs)
save_predicted_pct_misclassified_3panel()
save_predicted_pct_unclassified_3panel()
save_forecast_improved_misclassification(hybas_misclass_rate_maps$i/hybas_misclass_rate_maps10$i)
err_scatter$i | err_scatter$a | err_scatter$c
# ggsave("figures/iac_vs_pct_ghosts.png", width = 9, height = 3, dpi = 400)