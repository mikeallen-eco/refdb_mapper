# function to run leave-one-sequence-out model-based taxonomic assignment for a reference database 

# load libraries
library(dplyr)
library(tidyr)
library(Biostrings)

LOSO_models <- function(refdb,
                          out,
                          min_n = 2,
                          start_seq = 1,
                          random_seed = 100,
                          verbose = F){
  # set paths to RDP reference database for testing
  testmode <- F
  if(testmode %in% T){
    refdb <- refdb_Vences_16S
    out <- paste0(dirname(refdb_Vences_16S),"/loso_rdp/")
    min_n = 2
    random_seed = 100
    verbose = T
    start_seq = 1
  }
  
  # hard coded settings
  minBoot <- 0
  
  # read in arguments to pass through to helper functions
  RDP_out <- out
  
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
  
  # dataframe of reference db
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  # get vector of sequences with n >= min_n to loop through in broader database
  seqnums_vector <- which(rn$n >= min_n)
  seq_range <- start_seq:length(seqnums_vector)
  
  # loop through species list to test for
  for (i in start_seq:length(test_seqnums_vector)) { 
    message("Testing sequence ", i, " of ", length(test_seqnums_vector), ": ", rn[seqnums_vector[i],]$s)
    
    seq_num <- seqnums_vector[i]
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOSO_subset(refdb = r, seq_num, return_db = TRUE)
    
    # create a reference database of just the LOO target sequence
    refdb_looseq <- LOSO_subset(refdb = r, seq_num, return_db = FALSE)
    
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
              file = paste0(out, "RDP_LOSO_", 
                            RDP_df$qseqid[1],".csv"),
              row.names = F)
    
    # BA_df <- run_BayesANT(q_seqs = seqs,
    #                       ref_seqs = refdb_dim,
    #                       out = BA_out) %>%
    #   mutate(tmp = qseqid) %>%
    #   separate(tmp, into = c("gen", "sp", "acc"), sep = "_") %>%
    #   mutate(target_sp = paste0(gen, "_", sp)) %>%
    #   select(-gen, -sp)
    # 
    # write.csv(BA_df, 
    #           file = paste0(out, "BA_LOOseq_data_", 
    #                         RDP_df$qseqid[1],".csv"),
    #           row.names = F)
  }
  
}