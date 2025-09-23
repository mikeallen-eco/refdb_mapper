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

# Build synonyms list for phylogeny names
# phyl_synonyms <- build_rgbif_synonym_database(phyl_harmonized$uid)
# write.csv(phyl_synonyms, "data/phylogeny_rgbif_synonyms.csv", row.names = F)
phyl_synonyms <- read.csv("data/phylogeny_rgbif_synonyms.csv") %>%
  mutate(Accepted = ununderscore(Accepted), Synonym = ununderscore(Synonym))

# MOL to phylogeny
mol_format <- mol %>% select(uid = Accepted) %>% distinct() %>% mutate(uid = underscore(uid),
                                                                       full_sci_name = ununderscore(uid),
                                                                       genus_species = ununderscore(uid))
mol_to_phyl_harmonized <- harmonize_with_backbone(query = mol_format,
                                                  backbone = phyl_synonyms,
                                                  fuzzy_threshold = 0.97,
                                                  manual_tax = "data/mol_to_phyl_mammals_manual_notes.tsv")
write.csv(mol_to_phyl_harmonized, "data/mol_to_phyl_mammals_harmonized.csv", row.names = F)

# ---- Step 2 tally seqs per species & marker & join to hydrobasins

hydrobasin_refdb_info <- tally_sequences(hydrobasin_species_path, 
                                            refdb_cur_path = refdb_cur_paths, 
                                            refdb_harmonized_path)

# find widespread species that still need manual taxonomic matching (to prioritize)
top_n <- head(sort(table(hydrobasin_refdb_info$mol_name), decreasing = T), n = 500) %>% c()
to_go <- mol_to_phyl_harmonized[is.na(mol_to_phyl_harmonized$BB_Accepted),]$uid
top_n_and_to_go <- names(top_n) %in% to_go
(to_go_top_n <- to_go[top_n_and_to_go])

# ---- Step 3 get NND for each sp within hydrobasins & join to seq info
  
# hybas_nnd <- get_NDD_per_sp_all_hydrobasins(hydrobasin_map = hydrobasin_map, 
#                                             hydrobasin_species = hydrobasin_species_path) # ~ .74 s per hydrobasin
# saveRDS(hybas_nnd, "~/Documents/mikedata/refdb_mapper/hybas_nnd_world_20250922.rds")
hybas_nnd <- readRDS("~/Documents/mikedata/refdb_mapper/hybas_nnd_world_20250922.rds") 
hybas_nnd <- hybas_nnd %>%
  do.call(bind_rows, .) %>%
  select(HYBAS_ID, mol_name, phyl_name, nnd)

hydrobasin_refdb_nnd_info <- hydrobasin_refdb_info %>%
  left_join(hybas_nnd,
            by = join_by(HYBAS_ID, mol_name)) %>%
  select(HYBAS_ID, order, family, mol_name, phyl_name, nnd, contains("12S"), contains("16S")) %>%
  arrange(HYBAS_ID, order, nnd)
rm(hydrobasin_refdb_info)

hydrobasin_ref_nnd_info_sum <- summarize_by_hydrobasin(hydrobasin_refdb_nnd_info)

# ---- Step 4 map all summary stats by hydrobasin

# map results
to_plot = c("pct_ghosts", "med_seqs", "nnd_med", "total_species")
ghost_plot <- map_ghosts(df = hydrobasin_ref_nnd_info_sum,
                          hydrobasin_map = hydrobasin_map,
                          markers = markers,
                          to_plot = to_plot)
ghost_plot$Vences_16S$pct_ghosts
ghost_plot$MiMammalU_12S$pct_ghosts

# ---- Step 6 LOSO/LOSpO analysis (takes hours, saves csv files for convenience)

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
                  out = paste0(dirname(refdb_Mamm01_12S),"/"), start_seq = 1028)

LOSpO_ghostblaster(refdb = refdb_Mamm01_12S,
                   out = paste0(dirname(refdb_Mamm01_12S),"/"), start_seq = 1)

# Taylor_16S (92 bp)
LOSO_ghostblaster(refdb = refdb_Taylor_16S,
                  out = paste0(dirname(refdb_Taylor_16S),"/"), start_seq = 1)

LOSpO_ghostblaster(refdb = refdb_Taylor_16S,
                   out = paste0(dirname(refdb_Taylor_16S),"/"), start_seq = 1)

# ---- Step 6 - get NND for each species WITHIN each refdb

NND_per_sp_within_refdb <- get_NND_per_sp_within_refdb(
  refdb_cur = refdb_cur_paths,
  refdb_harmonized = refdb_harmonized_path,
  mol_to_phyl_harmonized = mol_to_phyl_harmonized_path,
  extinct = ncbi_extinct
)

NND_per_sp_within_refdb <-
  lapply(1:length(markers), function(i) {
    NND_per_sp_within_refdb[[i]] %>%
      mutate(marker = markers[i])
  }) %>%
  do.call(bind_rows, .)

n_seqs_per_sp <- 
    hydrobasin_refdb_nnd_info %>%
      select(mol_name, phyl_name, contains("_16"), contains("_12S")) %>%
      distinct()

NND_per_sp_within_refdb_seqs <- 
  NND_per_sp_within_refdb %>%
  left_join(n_seqs_per_sp %>% select(-phyl_name),
            by = join_by(mol_name)) %>%
  mutate(across(.cols = c(contains("_12S"), contains("_16S")), 
                function(x) replace_na(list(x = 0))))

# FYI species that still had NA sequence totals (it seems because they had no sequences)...
# [1] "Alticola_macrotis"         "Bos_frontalis"             "Bos_grunniens"            
# [4] "Bos_indicus"               "Bos_primigenius"           "Bos_taurus"               
# [7] "Brotomys_voratus"          "Bubalus_bubalis"           "Caloprymnus_campestris"   
# [10] "Camelus_bactrianus"        "Camelus_dromedarius"       "Capra_hircus"             
# [13] "Equus_asinus"              "Equus_caballus"            "Felis_catus"              
# [16] "Geocapromys_thoracatus"    "Hippotragus_leucophaeus"   "Hydrodamalis_gigas"       
# [19] "Lama_glama"                "Lama_pacos"                "Lama_vicugna"             
# [22] "Leptonychotes_weddellii"   "Lobodon_carcinophaga"      "Neomonachus_schauinslandi"
# [25] "Neomonachus_tropicalis"    "Neotoma_floridana"         "Ommatophoca_rossii"       
# [28] "Otaria_flavescens"         "Ovis_aries"                "Palaeopropithecus_ingens" 
# [31] "Perameles_eremiana"        "Potorous_platyops"         "Pteropus_niger"           
# [34] "Pteropus_rodricensis"      "Pteropus_tokudae"          "Taeromys_dominator"       
# [37] "Thylacinus_cynocephalus"   "Tonatia_saurophila"        "Xenothrix_mcgregori"      
# [40] "Zalophus_japonicus"        "Melomys_rubicola"       

# ---- Step 7 - compile LOSO/LOSpO outcomes

get_loo_outcomes()


###### ---- plots

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
