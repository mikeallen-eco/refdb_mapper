# ---- SET-UP
source("R/setup.R")
source("R/settings_Vences_16S.R")

# ---- Step 0 Make taxon-specific reference database for LOSO analysis (takes ~ 10 min)

# curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
#                  out = out_path, l = l, L = L, db_name = db_name)

# --- Step 1 Harmonize phyogeny & ref db taxonomy to MOL

data_env <- prepare_target_data_for_harmonization(mol_tax = mol_tax,
                                      mol_group = "Mammals",
                                      tree_names = tree_names,
                                      refdb_cur = refdb_cur)
                                      
refdb_harmonized <- harmonize_with_mol(mol = data_env$mol,
                                      target = data_env$refdb,
                                      fuzzy_threshold = 0.97,
                                      manual_tax = manual_tax_refdb)

write.csv(refdb_harmonized, "data/refdb_harmonized.csv", row.names = F)

# phyl_harmonized <- harmonize_with_mol(mol = data_env$mol,
#                                       target = data_env$phyl,
#                                       fuzzy_threshold = 0.97,
#                                       manual_tax = NULL)

# load mol names for checking
mol <- data_env$mol

# ---- Step 2 identify & map ghosts

hydrobasin_ghosts <- identify_ghosts(hydrobasin_species, 
                                     refdb_cur_path = refdb_cur_path, 
                                     refdb_harmonized_path)

(hydrobasin_ghosts_tally <- tally_ghosts(hydrobasin_ghosts))

seq_info_summarized_by_hydrobasin <- summarize_by_hydrobasin(hydrobasin_ghosts)


# ---- Step 2 LOSO analysis (takes hours, saves csv files for convenience)

# LOSO_ghostblaster(LOSO_refdb = LOO_refdb_path,
#                   out = out_path, marker = marker, start_seq = 1)

# LOSpO_ghostblaster(LOSpO_refdb = LOO_refdb_path,
#                   out = out_path, marker = marker, start_seq = 1)



# ---- Step 3 fit model of error rates based on LOSO analysis, NND, and n_seqs (takes a few min)

NND_per_sp_within_refdb <- get_NND_per_sp_within_refdb(refdb = LOO_refdb_path)

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

# ---- Step 5 correlate ghosts & NDD

get_pct_ghost_residuals <- function(hybas_pred_map_sf = hybas_pred_map_sf){
mod_data <- hybas_pred_map_sf %>% 
  filter(!is.na(pct_ghosts),
         !is.na(mean_i))
mod <- lm(mean_i ~ pct_ghosts, data = mod_data)
plot(residuals(mod) ~ mod_data$mean_nnd)
plot(residuals(mod) ~ mod_data$med_seqs)

mod2 <- lm(mean_i ~ pct_ghosts + mean_nnd + med_seqs, data = mod_data)
summary(mod2)

rmap_data <- mod_data %>%
  mutate(r = residuals(mod))

return(rmap_data)
}

# ---- Step 10 make maps and plots

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
