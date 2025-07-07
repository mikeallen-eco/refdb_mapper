library(data.table)
library(ghostblaster)
library(dplyr)
library(Biostrings)
the_files <- list.files("R", full.names = T)
the_files_00 <- the_files[grepl(the_files, pattern = "00")]
the_files_help <- the_files[grepl(the_files, pattern = "help")]
lapply(the_files_00, FUN = source)
lapply(the_files_help, FUN = source)

# ---- SET-UP
refdb <- "~/Documents/mikedata/refdb_geo/MIDORI2_UNIQ_NUC_GB265_lrRNA_RDP.fasta"
taxon <- "Mammalia"
fwd <- "AGACGAGAAGACCCYdTGGAGCTT"
rev <- "GATCCAACATCGAGGTCGTAA"
out <- "~/Documents/mikedata/refdb_geo/working/"
dl_tax <- FALSE
db_name <- "V16S_mammalia_midori265_tax20250609"
l = 100
L = 300
start_seq = 1
marker = "Vences_16S"
min_n = 2
verbose = TRUE
BLAST_args = "-max_target_seqs 5000 -perc_identity 80 -qcov_hsp_perc 90"

# ---- Step 0 Make taxon-specific reference database for LOSO analysis
  # wrap these in a curate_amplicons function

taxon_refdb <- subset_raw_refdb_by_taxon(refdb, taxon) 

if(dl_tax == TRUE){
  download_NCBI_taxonomy_crabs(out = out)
}

extract_amplicons_crabs(ref_seqs = taxon_refdb, fwd = fwd, rev = rev, 
                                     out = out, verbose = TRUE)

clean_amplicons_crabs(input = "crabs_amplicons_pga.txt",
                                  out, l, L, db_name)

LOSO_refdb <- readDNAStringSet(paste0(out, "refdb_", db_name, ".fasta")) # 7745 seqs

# ---- Step 2 identify & map ghosts

hg <- identify_ghosts(
  hydrobasin_mammals_path = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv",
  geotax_path = "data/geotax.csv",
  LOSO_refdb_path = paste0(out, "refdb_", db_name, ".fasta")
)

(hg_tally <- tally_ghosts(hg))

ghost_map <- map_ghosts(df = hg,
                        hydrobasin_map = "~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg",
                        just_US = TRUE,
                        save_maps = TRUE)

# ---- Step 2 Make taxon-specific reference database for LOSO analysis

LOSO_ghostblaster(LOSO_refdb = paste0(out, "refdb_", db_name, ".fasta"),
                  out,
                  marker = "Vences_16S",
                  min_n = 2,
                  start_seq = 2003,
                  BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                  verbose = T)

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


