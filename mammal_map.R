library(data.table)
library(ghostblaster)
lapply(list.files("R", full.names = T), FUN = source)

rich <- fread("data/hybas_L6_mammal_genus_richness.csv")
all <- fread("data/hybas_L6_mammal_intersections_harmonized.csv")
length(unique(all$HYBAS_ID)) # 16381
length(unique(all$sciname)) # 6346
length(unique(all$order))

summary(rich$genus_richness)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00   29.00   42.00   51.06   67.00  168.00

## set up
refdb <- "~/Documents/mikedata/refdb_geo/MIDORI2_UNIQ_NUC_GB265_lrRNA_RDP.fasta"
taxon <- "Mammalia"
fwd <- "AGACGAGAAGACCCYdTGGAGCTT"
rev <- "GATCCAACATCGAGGTCGTAA"
out <- "~/Documents/mikedata/refdb_geo/working/"
dl_tax <- FALSE
db_name <- "V16S_midori265_20250609"
l = 100
L = 300

# ---- Step 1 Make taxon-specific reference database for LOSO analysis

taxon_refdb <- subset_raw_refdb_by_taxon(refdb, taxon) 

if(dl_tax == TRUE){
  download_NCBI_taxonomy_crabs(out = out)
}

# todo: split this into multiple functions
extract_amplicons_crabs(ref_seqs = taxon_refdb, fwd = fwd, rev = rev, 
                                     out = out)

clean_amplicons_crabs(input = "crabs_amplicons_pga.txt",
                                  out, l, L, db_name)

refdb <- readDNAStringSet(paste0(out, "refdb_V16S_midori265_20250609.fasta")) # 7745 seqs

# ---- Step 2 Make taxon-specific reference database for LOSO analysis

prepare_data_for_LOSO()

# loop through species list to test for
# allsp_results_list <- list()
for (i in start_seq:length(local_seqnums_vector)) { 
  message("Testing sequence ", i, " of ", length(local_seqnums_vector), ": ", rn[local_seqnums_vector[i],]$s)
  tictoc::tic()
  
  seq_num <- local_seqnums_vector[i]
  
  # create diminished reference database subset (minus target sequence)
  refdb_dim <- LOOseq_refdb_subset(refdb, seq_num, return_db = TRUE)
  
  # write diminished reference database as a fasta in a tmp folder
  writeXStringSet(refdb_dim, 
                  filepath = paste0(GB_out,"tmp_GB_refdb.fasta"),
                  append = F)
  
  # create a reference database of just the LOO target sequence
  refdb_looseq <- LOOseq_refdb_subset(refdb, seq_num, return_db = FALSE) # , verbose = T
  
  # format the target sequence for taxonomy assignment
  seqs <- as.character(refdb_looseq)
  
  # run ghostblaster
  ghost_data <- ghostblaster(seqs,
                             refdb = paste0(GB_out, "tmp_GB_refdb.fasta"),
                             out = GB_out,
                             locals = locals,
                             ectoPredict = ectoPredict,
                             verbose = verbose,
                             BLAST_args = BLAST_args) 
  
  write.csv(ghost_data, 
            file = paste0(GB_out, "ghost_data_", 
                          ghost_data$qseqid[1],".csv"),
            row.names = F)
  
  tictoc::toc()
  
}

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
  # % of species with ≥ n sequences
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


