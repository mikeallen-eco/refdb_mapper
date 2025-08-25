# function to run leave-one-species-out loop for ghostblaster

LOSpO_ghostblaster <- function(LOSpO_refdb,
                                out,
                                marker,
                                start_seq = 1,
                                blast_path = "/usr/local/ncbi/blast/bin",
                                BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                                verbose = TRUE){

  # create directories if needed
  if (!dir.exists(out)) {
    dir.create(out)
  }
  
  # make a temporary directory for diminished reference databases
  if (!dir.exists(paste0(out, "lospo/"))) {
    dir.create(paste0(out, "lospo/"))
  }
  
  # read in reference database
  if (grepl(LOSpO_refdb[1], pattern = "fasta")) {
    r <- readDNAStringSet(LOSpO_refdb)
  } else{
    r <- LOSpO_refdb
  }
  
  # get df of reference database
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  sp_names <- sort(unique(rn$s))
  
  # make a locals list placeholder (nothing is local) until locals = NULL is implemented
  loc_p <- c("Yeti yeti", "Gallus gallus")
  
  # loop through species list to test for
  # allsp_results_list <- list()
  for (i in start_seq:length(sp_names)) { 
    if(verbose == TRUE){message("Testing sequence ", i, " of ", length(sp_names), ": ", sp_names[i])}
    tictoc::tic()
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOSpO_subset(refdb = LOSpO_refdb,
                              species = sp_names[i], 
                              return_db = TRUE)
    
    # isolate the LOSpO target sequence
    refdb_loosp <- LOSpO_subset(refdb = LOSpO_refdb,
                                species = sp_names[i], 
                                return_db = FALSE)
    
    # run ghostblaster
    ghost_data <- ghostblaster(as.character(refdb_loosp),
                               refdb = refdb_dim,
                               out = out,
                               locals = loc_p,
                               marker = marker,
                               verbose = F,
                               BLAST_args = BLAST_args) 
    
    write.csv(ghost_data$ghost_data, 
              file = paste0(out, "lospo/gd_lospo_", i, "_",
                            ghost_data$ghost_data$qseqid[1],".csv"),
              row.names = F)
    
    tictoc::toc()
    
  }
  
}