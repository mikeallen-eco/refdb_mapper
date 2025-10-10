# function to run leave-one-sequence-out model-based taxonomic assignment for a reference database 

# load libraries
library(dplyr)
library(tidyr)
library(Biostrings)

LOSO_models <- function(refdb,
                          out,
                          min_n = 2,
                          start_seq = 1,
                          extinct = ncbi_extinct,
                          random_seed = 100,
                          verbose = F){
  
  # read in arguments to pass through to helper functions
  RDP_out <- out
  
  # create directories if needed
  if (!dir.exists(file.path(RDP_out, "loso_rdp"))) {
    dir.create(file.path(RDP_out, "loso_rdp"))
  }
  
  # read in broader reference database fasta or DNA stringset
  if(grepl(refdb, pattern = "fasta")[1]){
    r <- readDNAStringSet(refdb)
  }else{r <- refdb}
  
  # Extract species names
  nms <- names(r)
  species <- sub(".*;", "", nms)
  keep <- !(species %in% extinct | grepl("_x_", species))
  
  # Subset the DNAStringSet to exclude extinct & hybrid species
  r_filtered <- r[keep]
  
  # get df of broader reference database
  rn <- RDP_to_dataframe(r_filtered)
  
  if(verbose %in% T){
    check_taxonomy_consistency(rn)
  }
  
  # dataframe of reference db
  rn <- rn %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  # get vector of sequences with n >= min_n to loop through in broader database
  seqnums_vector <- which(rn$n >= min_n)
  seq_range <- start_seq:length(seqnums_vector)
  
  # loop through species list to test for
  for (i in start_seq:length(seqnums_vector)) { 
    message("Testing sequence ", i, " of ", length(seqnums_vector), ": ", rn[seqnums_vector[i],]$s)
    
    seq_num <- seqnums_vector[i]
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOSO_subset(refdb = r_filtered, seq_num, return_db = TRUE)
    
    # create a reference database of just the LOO target sequence
    refdb_looseq <- LOSO_subset(refdb = r_filtered, seq_num, return_db = FALSE)
    
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
              file = file.path(out, "loso_rdp", paste0("RDP_loso_", i, "_", RDP_df$qseqid[1], ".csv")),
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