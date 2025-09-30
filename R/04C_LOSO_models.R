# function to run leave-one-sequence-out model-based taxonomic assignment for a reference database 

# load libraries
library(dplyr)
library(tidyr)
library(Biostrings)
source("R/functions/LOOseq_refdb_subset.R")
# source("R/functions/loo_RDP_test_single_seq.R")
source("R/functions/subset_refdb.R")
source("R/functions/RDP_to_dataframe.R")
source("R/functions/RDP_to_dada2.R")
source("R/functions/RDP_to_BayesANT.R")
source("R/functions/run_RDP.R")
source("R/functions/run_BayesANT.R")

LOOseq_models <- function(refdb, # a reference database including the locals and more (e.g., all in genera or all regional species)
                          locals, # set of 'local' species to include in testing (should rename this 'mock_community')
                          out,
                          start_seq = 1,
                          random_seed = 100,
                          verbose = F){
  # set paths to RDP reference database for testing
  testmode <- F
  if(testmode %in% T){
    refdb <- "data/midori_GB263_12SV5.pga70.uni.l85.L130.n10.rdp.ama.fasta"
    # missing families in amagen:
    # Oceanites, Myzornis, Aleadryas, Colluricincla, Daphoenositta, Eulacestoma, Falcunculus, Pachycephala, Rhagologus, Turnagra, Platysteira 
    locals <- read.csv("data/sp_names_tax.CdT20240712.csv")$orig_name
    # test_broader_db_only <- T
    out <- "data/LOOseq_models_12SV5_data/"
    # method = "rdp"
    random_seed = 100
    verbose = T
    library(dplyr)
    library(tidyr)
    library(Biostrings)
    source("R/functions/LOOseq_refdb_subset.R")
    # source("R/functions/loo_RDP_test_single_seq.R")
    source("R/functions/subset_refdb.R")
    source("R/functions/RDP_to_dataframe.R")
    source("R/functions/RDP_to_dada2.R")
    source("R/functions/RDP_to_BayesANT.R")
    source("R/functions/run_RDP.R")
    source("R/functions/run_BayesANT.R")
  }
  
  # hard coded settings
  minBoot <- 0
  
  # read in arguments to pass through to helper functions
  # needed?
  local_MOL_names <- locals
  RDP_out <- BA_out <- out
  
  # create directories if needed
  if (!dir.exists(RDP_out)) {
    dir.create(RDP_out)
  }
  
  # read in broader reference database fasta or DNA stringset
  if(grepl(refdb, pattern = "fasta")[1]){
    r <- readDNAStringSet(refdb)
  }else{r <- refdb}
  
  # get df of broader reference database
  rn <- RDP_to_dataframe(r)
  
  if(verbose %in% T){
    check_taxonomy_consistency(rn)
  }
  
  # get refdb subset to Tumbira species 
  r.loc <- subset_refdb(RDPstringset = r,
                        locals = local_MOL_names)
  
  # get Tumbira species names
  rn.loc <- RDP_to_dataframe(r.loc) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup() %>%
    filter(n > 1)
  sp_names.loc <- unique(rn.loc$s)
  
  # get vector of local sequences (with n > 1) to loop through in broader database
  local_seqnums_vector <- grep(paste(sp_names.loc, 
                                     collapse = "|"), names(r))
  
  # loop through species list to test for
  # allsp_results_list <- list()
  for (i in start_seq:length(local_seqnums_vector)) { 
    message("Testing sequence ", i, " of ", length(local_seqnums_vector), ": ", rn[local_seqnums_vector[i],]$s)
    
    seq_num <- local_seqnums_vector[i]
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOOseq_refdb_subset(refdb, seq_num, return_db = TRUE)
    
    # create a reference database of just the LOO target sequence
    refdb_looseq <- LOOseq_refdb_subset(refdb, seq_num, return_db = FALSE) # , verbose = T
    
    # format the target sequence for taxonomy assignment
    seqs <- as.character(refdb_looseq)
    
    RDP_df <- run_RDP(q_seqs = seqs,
                      ref_seqs = refdb_dim,
                      out = RDP_out,
                      seed = random_seed) %>%
      mutate(tmp = qseqid) %>%
      separate(tmp, into = c("gen", "sp", "acc"), sep = "_") %>%
      mutate(target_sp = paste0(gen, "_", sp)) %>%
      select(-gen, -sp)
    
    write.csv(RDP_df, 
              file = paste0(out, "RDP_LOOseq_data_", 
                            RDP_df$qseqid[1],".csv"),
              row.names = F)
    
    BA_df <- run_BayesANT(q_seqs = seqs,
                          ref_seqs = refdb_dim,
                          out = BA_out) %>%
      mutate(tmp = qseqid) %>%
      separate(tmp, into = c("gen", "sp", "acc"), sep = "_") %>%
      mutate(target_sp = paste0(gen, "_", sp)) %>%
      select(-gen, -sp)
    
    write.csv(BA_df, 
              file = paste0(out, "BA_LOOseq_data_", 
                            RDP_df$qseqid[1],".csv"),
              row.names = F)
  }
  
}