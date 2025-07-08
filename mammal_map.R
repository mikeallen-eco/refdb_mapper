library(data.table)
library(ghostblaster)
library(dplyr)
library(Biostrings)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00|01|02|help")]
lapply(the_files_to_load, FUN = source)

# ---- SET-UP
out_path <- "~/Documents/mikedata/refdb_geo/mammals_Vences_16S/"
db_name <- "V16S_mammalia_midori265_tax20250609"

# ---- Step 0 Make taxon-specific reference database for LOSO analysis

curate_amplicons(refdb = "~/Documents/mikedata/refdb_geo/MIDORI2_UNIQ_NUC_GB265_lrRNA_RDP.fasta", 
                 taxon = "Mammalia", 
                 fwd = "AGACGAGAAGACCCYdTGGAGCTT", 
                 rev = "GATCCAACATCGAGGTCGTAA",
                 out = out_path, 
                 l = 100, L = 300, 
                 db_name = db_name,
                 dl_tax = FALSE,
                 crabs_dir = "/Users/mikea/miniconda3/bin/conda", 
                 crabs_env = "crb2", 
                 verbose = TRUE)

# ---- Step 1 identify & map ghosts

hg <- identify_ghosts(
  hydrobasin_species = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv",
  geotax_path = "data/geotax.csv",
  LOSO_refdb_path = paste0(out_path, "refdb_", db_name, ".fasta")
)

(hg_tally <- tally_ghosts(hg))

ghost_map <- map_ghosts(df = hg,
                        hydrobasin_map = "~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg",
                        just_US = TRUE,
                        save_maps = TRUE)

# ---- Step 2 LOSO analysis

LOSO_ghostblaster(LOSO_refdb = paste0(out_path, "refdb_", db_name, ".fasta"),
                  out = out_path,
                  marker = "Vences_16S",
                  min_n = 2,
                  start_seq = 1,
                  BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                  verbose = TRUE)

# ---- Step 3 fit model of error rates based on LOSO analysis, NND, and order

loso_refdb_ndd <- get_NND_per_sp_within_refdb(
  refdb = paste0(out_path, "refdb_", db_name, ".fasta"),
  min_n = 2,
  # only calculate for species with at least min_n sequences
  tree = "data/phyl.tre",
  phyltax = "data/phyltax.csv",
  verbose = TRUE
)

loso_GB_compiled <- compile_loso_ghostblaster_results(loso_path = paste0(out_path, "loso"),
                           loso_refdb_ndd_df = loso_refdb_ndd)

loso_GB_outcomes <- compile_loso_ghostblaster_outcomes(loso_GB_compiled)


# model_blast_error_rates <- function(loso_GB_outcomes){

library(brms) 
blast_df <- loso_GB_outcomes$b %>%
  group_by(true_ncbi_name, order, nnd) %>%
  summarize(thresh98 = sum(thresh98_i),
            ecotag_i = sum(ecotag_i),
            conf90_i = sum(conf90_i),
            thresh98_i = sum(thresh98_i),
            ecotag_c = sum(ecotag_c),
            thresh99_i = sum(thresh99_i),
            thresh98_c = sum(thresh98_c),
            conf90_c = sum(conf90_c),
            n = length(true_ncbi_name),
            .groups = "drop")

library(lme4)

m_ecotag_i <- glmer(
  cbind(ecotag_i, n - ecotag_i) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_ecotag_i)

m_thresh98_i <- glmer(
  cbind(thresh98_i, n - thresh98_i) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_thresh98_i)

m_thresh99_i <- glmer(
  cbind(thresh99_i, n - thresh98_i) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_thresh99_i)

m_conf90_i <- glmer(
  cbind(conf90_i, n - conf90_i) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_conf90_i)

m_ecotag_c <- glmer(
  cbind(ecotag_c, n - ecotag_c) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_ecotag_c)

m_thresh98_c <- glmer(
  cbind(thresh98_c, n - thresh98_c) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_thresh98_c)

m_conf90_c <- glmer(
  cbind(conf90_c, n - conf90_i) ~ nnd + (1 | true_ncbi_name),
  family = binomial,
  data = blast_df
); summary(m_conf90_c)

true_ncbi_name <- unique(loso_GB_outcomes$b$true_ncbi_name)

preds_ecotag_i <- predict(
  m_ecotag_i,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_thresh98_i <- predict(
  m_thresh98_i,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_thresh99_i <- predict(
  m_thresh99_i,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_conf90_i <- predict(
  m_conf90_i,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
    # note: do not need to specify true_ncbi_name here if excluding random effects
  ),
  type = "response",
  re.form = NA
)

preds_ecotag_c <- predict(
  m_ecotag_c,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
  ),
  type = "response",
  re.form = NA
)


preds_thresh98_c <- predict(
  m_thresh98_c,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_conf90_c <- predict(
  m_conf90_c,
  newdata = expand.grid(
    nnd = seq(0, 200, length.out = 100)
    # note: do not need to specify true_ncbi_name here if excluding random effects
  ),
  type = "response",
  re.form = NA
)

pred_df_ecotag_i <- expand.grid(nnd = seq(0, 200, length.out = 100),
                       n = 100) %>%
  mutate(preds = preds_ecotag_i)

pred_df_thresh98_i <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                n = 100) %>%
  mutate(preds = preds_thresh98_i)

pred_df_thresh99_i <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                  n = 100) %>%
  mutate(preds = preds_thresh99_i)

pred_df_conf90_i <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                n = 100) %>%
  mutate(preds = preds_conf90_i)

pred_df_ecotag_c <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                n = 100) %>%
  mutate(preds = preds_ecotag_c)

pred_df_thresh98_c <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                  n = 100) %>%
  mutate(preds = preds_thresh98_c)

pred_df_conf90_c <- expand.grid(nnd = seq(0, 200, length.out = 100),
                                n = 100) %>%
  mutate(preds = preds_conf90_c)

ggplot(pred_df_ecotag_i) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,10))

ggplot(pred_df_thresh98_i) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,10))

ggplot(pred_df_thresh99_i) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,10))

ggplot(pred_df_conf90_i) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,10))

ggplot(pred_df_thresh98_c) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,100))

ggplot(pred_df_conf90_c) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,100))

ggplot(pred_df_ecotag_c) +
  geom_line(aes(x = nnd, y = 100*preds)) +
  scale_y_continuous(limits = c(0,100))

# FORMAT MODEL PREDICTION DATA FOR PLOTTING

# format plot data - LOOsp FP
set.seed(523); plotdata <- 
  tidybayes::epred_draws(b_ecotag_i, newdata = expand.grid(
    nnd = seq(0,200,length.out = 100),
    n = 100,
    order = orders), 
    probs = c(0.025, 0.1, 0.5, 0.90, 0.975)) %>%
  group_by(nnd, order) %>%
  summarize(mean = mean(.epred),
            Q2.5 = quantile(.epred, 0.025),
            Q10 = quantile(.epred, 0.1),
            Q50 = quantile(.epred, 0.5),
            Q90 = quantile(.epred, 0.9),
            Q97.5 = quantile(.epred, 0.975),
            .groups = "drop") %>% 
  mutate(metric = "incorrect",
         marker = "V16S",
         method = "BLAST (ecotag)")

library(ggplot2)
library(wesanderson)
col5 <- c(wes_palettes$FantasticFox1[c(1,2,3)], "darkgray", wes_palettes$Zissou1[5])

(plot <-
    plotdata %>%
    ggplot() +
    geom_ribbon(aes(x = nnd, 
                    ymin = Q2.5, ymax = Q97.5, 
                    fill = order),
                alpha = 0.25) + # color = "darkgray", 
    geom_line(aes(x = nnd, y = Q50,
                  color = order),
              linewidth = 1) +
    # facet_wrap(~method) +
    # facet_grid(cols = vars(method), rows = vars(type),
    #            scale = "free_y") +
    scale_x_continuous(limits = c(0,200)) +
    # scale_color_manual(values = col5) +
    # scale_fill_manual(values = col5) +
    labs(x = "Nearest evolutionary neighbor (MY)",
         y = "Probability incorrect (%)",
         color = "Taxonomic\ngroup") +
    theme_bw() +
    theme(text = element_text(size = 15),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.position = "bottom") +
    guides(fill = "none")
)

ggsave("figures/Fig4_plot_nearest_evgap_vs_LOO_accuracy_w_BLAST8.png",
       dpi = 400, width = 8, height = 5)

# ---- Step 4 get NND for each sp within hydrobasins for error rate predictions
  
  
  species_NND <- get_NND_per_sp_within_list(species_list)

# get nearest neighbor distance within watersheds for each species on list


# Map: % ghosts; or % of species with ≥ n (0,1,2,10) sequences
# mean n sequences per species
# model P(misclass) | NND & order for a novel sequence
# % of species likely to be misclassified (≥ X% probability); 
# avg. predicted misclassification rate etc.; 
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


