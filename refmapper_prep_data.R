# ---- SET-UP
source("R/setup.R")

# ---- Step 0 Make taxon-specific reference databases for LOSO analysis (takes ~ 10 min)

source("R/settings_Vences_16S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_Taylor_16S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_V5_12S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_MiMamm_12S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

source("R/settings_Mamm01_12S.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev,
                 out = out_path, l = l, L = L, db_name = db_name)

# --- Step 1 Harmonize phyogeny & ref db taxonomy to MOL

# ref db to MOL
refdb_harmonized <- harmonize_with_backbone(backbone = data_env$mol,
                                            query = data_env$refdb,
                                            fuzzy_threshold = 0.97,
                                            manual_tax = manual_tax_refdb)

write.csv(refdb_harmonized, "data/refdb_mammals_harmonized.csv", row.names = F)

# phylogeny to MOL
phyl_harmonized <- harmonize_with_backbone(backbone = data_env$mol,
                                           query = data_env$phyl,
                                           fuzzy_threshold = 0.97,
                                           manual_tax = manual_tax_phyl)

write.csv(phyl_harmonized, "data/phyl_mammals_harmonized.csv", row.names = F)

# MOL to phylogeny
# Build synonyms list for phylogeny names
# phyl_synonyms <- build_rgbif_synonym_database(phyl_harmonized$uid)
# write.csv(phyl_synonyms, "data/phylogeny_rgbif_synonyms.csv", row.names = F)
phyl_synonyms <- read.csv("data/phylogeny_rgbif_synonyms.csv") %>%
  mutate(Accepted = ununderscore(Accepted), Synonym = ununderscore(Synonym))

mol_format <- mol %>% select(uid = Accepted) %>% distinct() %>% mutate(uid = underscore(uid),
                                                                       full_sci_name = ununderscore(uid),
                                                                       genus_species = ununderscore(uid))
mol_to_phyl_harmonized <- harmonize_with_backbone(query = mol_format,
                                                  backbone = phyl_synonyms,
                                                  fuzzy_threshold = 0.94, # only a few wrong fuzzy matches to manually correct
                                                  manual_tax = "data/mol_to_phyl_mammals_manual_notes.tsv")
write.csv(mol_to_phyl_harmonized, "data/mol_to_phyl_mammals_harmonized.csv", row.names = F)

# ---- Step 2 tally seqs per species & marker & join to hydrobasins

hydrobasin_refdb_info <- tally_sequences(hydrobasin_species_path, 
                                            refdb_cur_path = refdb_cur_paths, 
                                            refdb_harmonized_path)

# find widespread species that still need manual taxonomic matching (to prioritize)
priority_species_to_harmonize(harmonized_df = mol_to_phyl_harmonized, show = 10)
priority_species_to_harmonize(harmonized_df = phyl_harmonized, show = 10)

# ---- Step 3 get NND for each sp within hydrobasins & join to seq info
  
hybas_nnd <- get_NDD_per_sp_all_hydrobasins_and_markers(hydrobasin_refdb_info,
                                                        markers = markers, n_cores = 4)

# hybas_nnd <- get_NDD_per_sp_all_hydrobasins(hydrobasin_refdb_info = hydrobasin_refdb_info[1:1000,],
#                                             markers = markers) # ~ .7 s per hydrobasin
saveRDS(hybas_nnd, "~/Documents/mikedata/refdb_mapper/hybas_nnd_world_all_markers_20250930.rds")
hybas_nnd <- readRDS("~/Documents/mikedata/refdb_mapper/hybas_nnd_world_all_markers_20250930.rds") 

# ---- Step 4 LOSO/LOSpO analysis (takes hours, saves csv files for convenience)

# Vences_16S (237 bp)
LOSO_ghostblaster(refdb = refdb_Vences_16S,
                  out = paste0(dirname(refdb_Vences_16S),"/"), start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_Vences_16S,
                   out = paste0(dirname(refdb_Vences_16S),"/"), start_seq = 1)

# RiazVert1_12S (105 bp)
LOSO_ghostblaster(refdb = refdb_RiazVert1_12S,
                  out = paste0(dirname(refdb_RiazVert1_12S),"/"), start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_RiazVert1_12S,
                  out = paste0(dirname(refdb_RiazVert1_12S),"/"), start_seq = 1)

# MiMammalU_12S (169 bp)
LOSO_ghostblaster(refdb = refdb_MiMammalU_12S,
                  out = paste0(dirname(refdb_MiMammalU_12S),"/"), start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_MiMammalU_12S,
                   out = paste0(dirname(refdb_MiMammalU_12S),"/"), start_seq = 1)

# Mamm01_12S (59 bp)
LOSO_ghostblaster(refdb = refdb_Mamm01_12S,
                  out = paste0(dirname(refdb_Mamm01_12S),"/"), min_length = 40, start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_Mamm01_12S,
                   out = paste0(dirname(refdb_Mamm01_12S),"/"), min_length = 40, start_seq = 1)

# Taylor_16S (92 bp)
LOSO_ghostblaster(refdb = refdb_Taylor_16S,
                  out = paste0(dirname(refdb_Taylor_16S),"/"), min_length = 40, start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_Taylor_16S,
                   out = paste0(dirname(refdb_Taylor_16S),"/"), min_length = 40, start_seq = 1)

# ---- Step 5 - build predictor data for error model

# get NND and n seqs for each species within each refdb
error_model_predictor_data <- build_error_model_data(
  refdb_cur = refdb_cur_paths,
  refdb_harmonized = refdb_harmonized_path,
  mol_to_phyl_harmonized = mol_to_phyl_harmonized_path,
  extinct = c(ncbi_extinct),
  marker_names = markers
)

# ---- Step 6 - compile outcomes for final error model data

mdat <- get_loo_outcomes(marker_directories = dirname(refdb_cur_paths)[1:5],
                 markers = markers[1:5],
                 refdb_harmonized = refdb_harmonized_path,
                 refdb_nnd = error_model_predictor_data)

# ---- Step 7 - fit models

fits <- fit_models_loso_lospo(assign_rubric = "thresh98",
                      markers = markers[1:5],
                      outcomes = mdat)

# ---- Step 8 - predict error rates for species within hydrobasins

complete_hybas_data <- predict_error_rate_hybas(hybas_info = hydrobasin_refdb_nnd_info,
                                                preds = fits,
                                                markers = markers[1:5])

# ---- Step 9 - make final map sf with data

final_sf <- make_complete_hybas_data_sf(complete_hybas_data = complete_hybas_data,
                                                      map = hydrobasin_map)

saveRDS(final_sf, "~/Documents/mikedata/refdb_mapper/final_hybas_data_sf_20250930.rds")
final_sf <- readRDS("~/Documents/mikedata/refdb_mapper/final_hybas_data_sf_20250930.rds")






# ... processes below here were (or will be) moved to refmapper_present_info.R script


# ---- Step 10 - return tables & info based on coordinates

hybas_info <- get_polygon_attributes_from_coords(hybas_data = final_sf)

# ideas
# marker combo that maximizes coverage
# marker combo that minimizes % incorrect 
# pairs of species most likely to be confused
# marker combo that maximizes % correct
  # for target sp list?


# ---- Step 10 - make model plots

eplots <- plot_predicted_loso_lopso_error(preds = fits, markers = markers[1:3])
eplots$RiazVert1_12S$loso$i
eplots$MiMammalU_12S$loso$i
eplots$Vences_16S$loso$i

# ---- Step 10 map all results by hydrobasin

# map results
to_plot = c("pct_ghosts", "med_seqs", "nnd_med", "total_species")
ghost_plot <- map_ghosts(df = hydrobasin_ref_nnd_info_sum,
                         hydrobasin_map = hydrobasin_map,
                         markers = markers,
                         to_plot = to_plot)
ghost_plot$Vences_16S$pct_ghosts
ghost_plot$MiMammalU_12S$pct_ghosts


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
