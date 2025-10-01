library(Biostrings)
library(dplyr)
library(tictoc)

LOSpO_models <- function(refdb,
                               out,
                               start_sp = 1,
                               random_seed = 100,
                               verbose = TRUE) {
  
  # read in arguments to pass through to helper functions
  RDP_out <- out
  
  # create directories if needed
  if (!dir.exists(RDP_out)) {
    dir.create(RDP_out)
  }
  
  if (!dir.exists(file.path(RDP_out, "lospo_rdp"))) {
    dir.create(file.path(RDP_out, "lospo_rdp"))
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
  
  # get df of reference database
  rn <- RDP_to_dataframe(r_filtered)
  sp_names <- sort(unique(rn$s))

  if(verbose %in% T){
    check_taxonomy_consistency(rn)
  }
  
  # loop through species list to test LOOsp for each
  tic()
  for (i in start_sp:length(sp_names)) { 
    species_name <- sp_names[i]
    message("Testing species ", i, " of ", length(sp_names), ": ", species_name)
    
    # diminished refdb (minus this species)
    refdb_dim <- LOSpO_subset(refdb = r_filtered, species = species_name, return_db = TRUE)
    # LOSpO target (just this species)
    refdb_loosp <- LOSpO_subset(refdb = r_filtered, species = species_name, return_db = FALSE)
    
    # format the target sequence for taxonomy assignment
    seqs <- as.character(refdb_loosp)
    
    RDP_df <- run_RDP(q_seqs = seqs,
                      ref_seqs = refdb_dim,
                      out = RDP_out,
                      seed = random_seed) %>%
      mutate(tmp = qseqid) %>%
      separate(tmp, into = c("gen", "sp", "acc"), sep = "_") %>%
      mutate(target_sp = paste0(gen, "_", sp)) %>%
      select(-gen, -sp)
    
    write.csv(RDP_df, 
              file = file.path(out, "lospo_rdp", paste0("RDP_lospo_", i, "_", species_name, ".csv")),
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
    #           file = paste0(out, "BA_LOOsp_data_", 
    #                         RDP_df$qseqid[1],".csv"),
    #           row.names = F)
    
  }

}
