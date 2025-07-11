library(data.table)
library(ghostblaster)
library(dplyr)
library(Biostrings)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00|01|02|03|04|05|help")]
lapply(the_files_to_load, FUN = source)

# ---- SET-UP
out_path <- "~/Documents/mikedata/refdb_geo/mammals_Vences_16S/"
db_name <- "V16S_mammalia_midori265_tax20250609"

# ---- Step 0 Make taxon-specific reference database for LOSO analysis
# note: add crabs subset command to subset final output to Mammalia (excludes random turtle & molusks that made it in)
curate_amplicons(refdb = "~/Documents/mikedata/refdb_geo/MIDORI2_UNIQ_NUC_GB265_lrRNA_RDP.fasta", 
                 taxon = "Mammalia", 
                 fwd = "AGACGAGAAGACCCYdTGGAGCTT", 
                 rev = "GATCCAACATCGAGGTCGTAA",
                 out = out_path, 
                 l = 100, L = 300, 
                 db_name = db_name,
                 dl_tax = FALSE,
                 conda_dir = "/Users/mikea/miniconda3/bin/conda", 
                 conda_env = "crb2", 
                 verbose = TRUE)

# ---- Step 1 identify & map ghosts

hydrobasin_ghosts <- identify_ghosts(
  hydrobasin_species = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv",
  geotax_path = "data/geotax.csv",
  LOSO_refdb_path = paste0(out_path, "refdb_", db_name, ".fasta")
)

(hydrobasin_ghosts_tally <- tally_ghosts(hydrobasin_ghosts))

hybas_subset <- subset_hydrobasins_countries(
  hydrobasin_map = "~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg",
  countries = c("United States of America", "Mexico", "Canada"))

(hydrobasin_ghosts_tally_sub <- tally_ghosts(hydrobasin_ghosts %>% filter(HYBAS_ID %in% hybas_subset$HYBAS_ID)))

(ghost_map <- map_ghosts(df = hydrobasin_ghosts,
                        hydrobasin_map = hybas_subset,
                        save_maps = "figures/",
                        suffix = "NA_no_AK"))

# ---- Step 2 LOSO analysis

LOSO_ghostblaster(LOSO_refdb = paste0(out_path, "refdb_", db_name, ".fasta"),
                  out = out_path,
                  marker = "Vences_16S",
                  min_n = 2,
                  start_seq = 1,
                  BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                  verbose = TRUE)

LOSpO_ghostblaster(LOSpO_refdb = paste0(out_path, "refdb_", db_name, ".fasta"),
                  out = out_path,
                  marker = "Vences_16S",
                  start_seq = 1,
                  BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                  verbose = TRUE)

# ---- Step 3 fit model of error rates based on LOSO analysis, NND, and n_seqs

NND_per_sp_within_refdb <- get_NND_per_sp_within_refdb(
  refdb = paste0(out_path, "refdb_", db_name, ".fasta"),
  min_n = 1,
  # only calculate for species with at least min_n sequences
  tree = "data/phyl.tre",
  phyltax = "data/phyltax.csv",
  verbose = TRUE
)

loo_outcomes <- get_loo_outcomes(loso_gb_path = paste0(out_path, "loso"),
                                  lopso_gb_path = paste0(out_path, "lospo"),
                                  refdb_nnd = NND_per_sp_within_refdb)

preds_loso_lospo <- get_preds_loso_lospo(assign = "BLAST",
                                         accuracy_metric = "thresh98",
                                         outcomes = loo_outcomes)

(predicted_loso_lopso_error_plot <- plot_predicted_loso_lopso_error(preds = preds_loso_lospo,
                                        save_to = "figures/",
                                        suffix = "BLAST98"))

# ---- Step 4 get NND, n_seqs, & error rate predictions for each sp within hydrobasins
  
hybas_subset <- subset_hydrobasins_countries(
  hydrobasin_map = "~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg",
  countries = c("United States of America", "Mexico", "Canada"))

# hybas_nnd <- get_NDD_per_sp_all_hydrobasins(
#   hydrobasin_map = hybas_subset,
#   hydrobasin_species = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv",
#   tree = "data/phyl.tre",
#   phyltax = "data/phyltax.csv",
#   sp_list_tax = "data/geotax.csv",
#   verbose = T
# ) # ~ .74 s per hydrobasin (1914 s for n = 2599 in NA)
# saveRDS(hybas_nnd, paste0(out_path, "hybas_nnd_NA.rds"))
hybas_nnd <- readRDS(paste0(out_path, "hybas_nnd_NA.rds"))

hybas_error_data <- format_hybas_error_data(hybas_ndd_df = hybas_nnd,
                                    ref_path = paste0(out_path,"refdb_", db_name, ".fasta"),
                                    preds_list = preds_loso_lospo)

# ---- Step 5 map error rates across hydrobasins

hybas_pred_map_sf <- make_hybas_pred_map_sf(hydrobasin_map = hybas_subset,
                                       hybas_pred = hybas_error_data)

(hybas_misclass_rate_maps <- plot_hybas_misclass_rate_maps(hybas_pred_map_sf = hybas_pred_map_sf, 
                                                          suffix = "BLAST98"))



# mean n sequences per species (non-ghosts)
# % of species likely to be misclassified (≥ X% probability); 

# exclude species from misclassification calcs that are likely unamplifiable
# generate lists of species w/ high likelihood of misclassifying or under-classifying etc.


# * e.g., model-based approaches, BLAST, or GhostBLASTer w/ various thresholds (I realized that it should be possible to calibrate probability scores for each hydroshed using a customized summarize function on the global leave-one-sequence-out data.)
# mean probability of misclass, underclass, correct for a novel sequence
# % of species with > X% P(correct), etc. 
# % of sequences misclassified in LOSO
# % of sequences underclassified in LOSO

# One LOSO analysis (using GhostBLASTer) for all mammals in ref db w/ > 1 seq (n_mam_seqs_gr1 ~ ?)
  # for (i in 1:n_mam_seqs_gr1)
  #   - leave one sequence out
  # - GhostBLAST against all others
  # - collect top hits (ghostdata) into database
  # - columns for qseqid & target_sp

# Perform LOSO species determination within each hydrobasin (n_hybas ~ 16000)
  # for (i in 1:n_hybas)
  # - subset LOSO ghostdata from mock community (union of ref spp & hybas i spp)
  # - filter to keep only hits for locals within hybas i
  # - generate ID and confidence metrics by qseqid for various BLAST rubrics, GhostBLASTer

# Map simple results
  # % of species with ≥ n (0,1,2,10) sequences
  # mean n sequences per species
  # % of sequences misclassified in LOSO
  # % of sequences underclassified in LOSO

# Use models to predict P(misclass), P(underclass), P(correct) of novel sequence
# collect all outcome data (incorrect / abstain / correct)
  # model P(misclass), P(underclass), P(correct) based on NND (MY)

# Map results based on modeled accuracy (takes ghosts into account)
  # mean probability of misclass, underclass, correct for a novel sequence
  # % of species with > X% P(correct), etc. 

# consider repeating for LNSO (leave-no-sequence-out)
  # interesting to see P(underclass) of species with exact variant present in database

# could also do for model based approaches (RDP, BayesANT)
  # but for now I can only imagine it being feasbile using GLOBAL referenece database for training
  # (otherwise you would have to train and test for each sequence / hydrobasin combo ~ > 80M iterations)


###### old

# Step 1 (later): make a better consensus tree for all mammals (all 10,000 trees instead of random 1000 trees)
# 
# Step 2: For each marker & taxonomic assignment rubric*
# 1. Get nearest evolutionary neighbor distance (NND in MY) for the set of mammal species that appear in MIDORI2 (considering only species in the database)
# 2. Do leave-one-sequence-out test for each mammal sequence in MIDORI2 (full database) - for each, collect real name, predicted name, confidence info [max pident, gap, conf, etc.]
# 3. Fit model to predict misclassification rate | order & NND
# 4. Reconcile Hybas6 mammal list taxonomy with NCBI and count available sequences for each
# 5. Reconcile Hybas6 mammal list taxonomy with phylogeny & get NND for each mammal within list (including ghosts) 
# 6. Use model to predict misclassification & under-classification rate for all non-ghosts | order & NND w/in watershed (including ghosts)
# 7. [optional: use a leave-one-species-out analysis to also predict misclassification rate for each ghost within each watershed]
# 8. Map: % ghosts; % of species likely to be misclassified (≥ X% probability); avg. predicted misclassification rate etc.; generate lists of species w/ high likelihood of misclassifying or under-classifying etc.
# * e.g., model-based approaches, BLAST, or GhostBLASTer w/ various thresholds (I realized that it should be possible to calibrate probability scores for each hydroshed using a customized summarize function on the global leave-one-sequence-out data.)


