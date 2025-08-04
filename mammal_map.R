# ---- SET-UP
source("R/setup.R")
source("R/settings_Vences_16S.R")

# ---- Step 0 Make taxon-specific reference database for LOSO analysis

curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, 
                 out = out_path, l = l, L = L, db_name = db_name)

# identify_nonamplifiers()

# ---- Step 1 identify & map ghosts

hydrobasin_ghosts <- identify_ghosts(LOO_refdb_path = LOO_refdb_path)

(hydrobasin_ghosts_tally <- tally_ghosts(hydrobasin_ghosts))

seq_info_summarized_by_hydrobasin <- summarize_by_hydrobasin(hydrobasin_ghosts)

(ghost_plots <- map_ghosts(df = seq_info_summarized_by_hydrobasin,
                        hydrobasin_map = hydrobasin_map))

(med_num_seqs_plot <- map_ghosts(df = seq_info_summarized_by_hydrobasin,
                           hydrobasin_map = hydrobasin_map, 
                           to_plot = "med_seqs"))

# ---- Step 2 LOSO analysis (takes hours, saves csv files for convenience)

LOSO_ghostblaster(LOSO_refdb = LOO_refdb_path,
                  out = out_path, marker = marker, start_seq = 1)

LOSpO_ghostblaster(LOSpO_refdb = LOO_refdb_path,
                  out = out_path, marker = marker, start_seq = 1)

# ---- Step 3 fit model of error rates based on LOSO analysis, NND, and n_seqs

NND_per_sp_within_refdb <- get_NND_per_sp_within_refdb(
  refdb = LOO_refdb_path)

loo_outcomes <- get_loo_outcomes(loso_gb_path = paste0(out_path, "loso"),
                                  lopso_gb_path = paste0(out_path, "lospo"),
                                  refdb_nnd = NND_per_sp_within_refdb)

preds_loso_lospo <- get_preds_loso_lospo(assign = "BLAST",
                                         accuracy_metric = "thresh98",
                                         outcomes = loo_outcomes)

predicted_loso_lospo_error_plots <- plot_predicted_loso_lopso_error(preds = preds_loso_lospo)

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

# ---- Step 5 map error rates & nnd across hydrobasins

hybas_misclass_rate_maps <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf)

hybas_mean_nnd_map <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf,
                                                    to_plot = "nnd")

# ---- Step 6 save plots to png files

save_num_mammals(ghost_plots$num_all)
save_mean_nnd(hybas_mean_nnd_map$nnd)
save_num_pct_molecular_ghosts(ghost_plots$num_ghosts/ghost_plots$pct_ghosts)
save_median_num_seqs_non_ghosts(med_seqs$med_seqs)
save_predicted_pct_misclassified_3panel()
save_predicted_pct_unclassified_3panel()
save_forecast_improved_misclassification(hybas_misclass_rate_maps$i/hybas_misclass_rate_maps10$i)

# ---- Step 7 correlate mean error rate with % ghosts & num spp

err_scatter <- error_rate_scatterplots(hybas_pred_map_sf)

err_scatter$i | err_scatter$a | err_scatter$c
# ggsave("figures/iac_vs_pct_ghosts.png", width = 9, height = 3, dpi = 400)

# ---- Step 8 example in one hydrobasin (Delaware River)

dewa <- -1529894232; rar <- -1529894602
# dewa_map <- hydrobasin_map %>% filter(HYBAS_ID %in% dewa)

# geotax <- fread("data/geotax.csv")
dewa_ghosts <- hydrobasin_ghosts %>%
  filter(HYBAS_ID %in% dewa)

(dewa_ghosts_tally <- tally_ghosts(dewa_ghosts))
dewa_info_summarized <- seq_info_summarized_by_hydrobasin %>%
  filter(HYBAS_ID %in% dewa)

dewa_hybas_error_data <- hybas_error_data %>%
  filter(HYBAS_ID %in% dewa)

dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

# ---- Step 9 map error rates across hydrobasins assuming ≥ 10 sequences per species

hybas_misclass_rate_maps10 <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf%>%
                                                              select(-starts_with("mean_")) %>%
                                                              rename_with(~ gsub("mean10", "mean", .),
                                                                          starts_with("mean10")))