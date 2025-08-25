# function to run leave-one-sequence-out loop for ghostblaster

LOSO_ghostblaster <- function(LOSO_refdb,
                                out,
                                marker,
                                min_n = 2,
                                start_seq = 1,
                                blast_path = "/usr/local/ncbi/blast/bin",
                                BLAST_args = "-max_target_seqs 5000 -perc_identity 75 -qcov_hsp_perc 90",
                                verbose = TRUE){

  # create directories if needed
  if (!dir.exists(out)) {
    dir.create(out)
  }
  
  # make a temporary directory for diminished reference databases
  if (!dir.exists(paste0(out, "loso/"))) {
    dir.create(paste0(out, "loso/"))
  }
  
  # read in reference database
  if (grepl(LOSO_refdb[1], pattern = "fasta")) {
    r <- readDNAStringSet(LOSO_refdb)
  } else{
    r <- LOSO_refdb
  }
  
  # get df of reference database
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  # get index positions of sequences from species with > min_n variants present
  seqnums_vector <- which(rn$n >= min_n)
  
  # make a locals list placeholder (nothing is local) until locals = NULL is implemented
  loc_p <- c("Yeti yeti", "Gallus gallus")
  
  # loop through species list to test for
  # allsp_results_list <- list()
  for (i in start_seq:length(seqnums_vector)) { 
    if(verbose == TRUE){message("Testing sequence ", i, " of ", length(seqnums_vector), ": ", rn[seqnums_vector[i],]$s)}
    tictoc::tic()
    
    seq_num <- seqnums_vector[i]
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOSO_subset(refdb = r, seq_num, return_db = TRUE)
    
    # isolate the LOSO target sequence
    refdb_looseq <- LOSO_subset(refdb = r, seq_num, return_db = FALSE)
    
    # run ghostblaster
    ghost_data <- ghostblaster(as.character(refdb_looseq),
                               refdb = refdb_dim,
                               out = out,
                               locals = loc_p,
                               marker = marker,
                               verbose = F,
                               BLAST_args = BLAST_args) 
    
    write.csv(ghost_data$ghost_data, 
              file = paste0(out, "loso/gd_", i, "_",
                            ghost_data$ghost_data$qseqid[1],".csv"),
              row.names = F)
    
    tictoc::toc()
    
  }
  
}