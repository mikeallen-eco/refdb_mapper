# ---- SET-UP
source("R/setup.R")

# ---- Step 0 Make taxon-specific reference databases for LOSO analysis (takes ~ 10 min)

source("R/settings_Vences_16S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_V5_12S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_MiMamm_12S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

# --- Step 1 Harmonize phyogeny & ref db taxonomy to MOL

refdb_cur <- c(refdb_V5_12S, refdb_MiMamm_12S, refdb_Vences_16S)
data_env <- prepare_target_data_for_harmonization(mol_tax = mol_tax,
                                      mol_group = "Mammals",
                                      tree_names = tree_names,
                                      refdb_cur = refdb_cur,
                                      extinct = extinct)
                                      
refdb_harmonized <- harmonize_with_mol(mol = data_env$mol,
                                      target = data_env$refdb,
                                      fuzzy_threshold = 0.97,
                                      manual_tax = manual_tax_refdb)

write.csv(refdb_harmonized, "data/refdb_mammals_harmonized.csv", row.names = F)

# phyl_harmonized <- harmonize_with_mol(mol = data_env$mol,
#                                       target = data_env$phyl,
#                                       fuzzy_threshold = 0.97,
#                                       manual_tax = NULL)

# load mol names for checking
mol <- data_env$mol

# ---- Step 2 identify & map ghosts

refdb_cur_paths <- c(V5_12S = refdb_V5_12S, MiMamm_12S = refdb_MiMamm_12S, Vences_16S = refdb_Vences_16S)
hydrobasin_refdb_info <- identify_ghosts(hydrobasin_species, 
                                            refdb_cur_path = refdb_cur_paths, 
                                            refdb_harmonized_path)

# ---- Step 4 get NND, n_seqs, & error rate predictions for each sp within hydrobasins
  
# hybas_nnd <- get_NDD_per_sp_all_hydrobasins() # ~ .74 s per hydrobasin
# saveRDS(hybas_nnd, paste0(out_path, "hybas_nnd_world_fixed.rds"))
hybas_nnd <- readRDS("~/Documents/mikedata/refdb_geo/hybas_nnd_world_fixed.rds") %>%
  do.call(bind_rows, .) %>%
  select(HYBAS_ID, sciname = geo_name, nnd)

hydrobasin_refdb_nnd_info <- hydrobasin_refdb_info %>%
  left_join(hybas_nnd,
            by = join_by(HYBAS_ID, sciname)) %>%
  select(HYBAS_ID, order, family, sciname, nnd, contains("12S"), contains("16S")) %>%
  arrange(HYBAS_ID, order, nnd)

# saveRDS(hydrobasin_refdb_nnd_info, "~/Documents/mikedata/refdb_geo/hydrobasin_refdb_nnd_info.rds")
hydrobasin_refdb_nnd_info <- readRDS("~/Documents/mikedata/refdb_geo/hydrobasin_refdb_nnd_info.rds")

hydrobasin_ref_nnd_info_sum <- summarize_by_hydrobasin(hydrobasin_refdb_nnd_info)

# ---- make maps and plots

(ghost_plot <- map_ghosts(df = seq_info_summarized_by_hydrobasin_MiMamm_12S,
                           hydrobasin_map = hydrobasin_map,
                           to_plot = "pct_ghosts"))

(med_num_seqs_plot <- map_ghosts(df = seq_info_summarized_by_hydrobasin,
                                 hydrobasin_map = hydrobasin_map,
                                 to_plot = "med_seqs"))

predicted_loso_lospo_error_plots <- plot_predicted_loso_lopso_error(preds = preds_loso_lospo)

hybas_misclass_rate_maps <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf)

hybas_mean_nnd_map <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf,
                                                    to_plot = "nnd")

err_scatter <- error_rate_scatterplots(hybas_pred_map_sf)

plot_pct_ghost_residuals <- function(rmap_data = rmap_data){
    
r_map <- rmap_data  %>%
    ggplot() +
    geom_sf(aes(fill = r),
            color = NA) +
    scale_fill_viridis_c(option = "inferno") +
    labs(fill = "%") +
    theme_minimal()

return(r_map)

}

pct_ghost_residuals_plot <- plot_pct_ghost_residuals(rmap_data)

# 
# plot forecasted improvement if all spp had ≥ 10 seqs
hybas_misclass_rate_maps10 <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf%>%
                                                              select(-starts_with("mean_")) %>%
                                                              rename_with(~ gsub("mean10", "mean", .),
                                                                          starts_with("mean10")))

# ---- Step 11 save maps and plots to png files

# save_num_mammals(ghost_plots$num_all)
# save_mean_nnd(hybas_mean_nnd_map$nnd)
# save_num_pct_molecular_ghosts(ghost_plots$num_ghosts/ghost_plots$pct_ghosts)
# save_median_num_seqs_non_ghosts(med_seqs$med_seqs)
# save_predicted_pct_misclassified_3panel()
# save_predicted_pct_unclassified_3panel()
# save_forecast_improved_misclassification(hybas_misclass_rate_maps$i/hybas_misclass_rate_maps10$i)
# err_scatter$i | err_scatter$a | err_scatter$c
# ggsave("figures/iac_vs_pct_ghosts.png", width = 9, height = 3, dpi = 400)
# save_pct_ghost_residuals(pct_ghost_residuals_plot = pct_ghost_residuals_plot)
