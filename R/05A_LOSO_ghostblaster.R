library(Biostrings)
library(dplyr)
library(tictoc)
library(parallel)

LOSO_ghostblaster <- function(refdb,
                              out,
                              marker = "average",
                              min_n = 2,
                              start_seq = 1,
                              blast_path = "/usr/local/ncbi/blast/bin",
                              BLAST_args = "-max_target_seqs 5000 -perc_identity 75",
                              min_length = 90,
                              verbose = TRUE,
                              parallel = FALSE,
                              n_cores = detectCores() - 1) {
  
  # create directories if needed
  if (!dir.exists(out)) dir.create(out)
  if (!dir.exists(file.path(out, "loso"))) dir.create(file.path(out, "loso"))
  if (!dir.exists(file.path(out, "tmp"))) dir.create(file.path(out, "tmp"))
  
  # read in reference database
  if (grepl(refdb[1], pattern = "fasta")) {
    r <- readDNAStringSet(refdb)
  } else {
    r <- refdb
  }
  
  # dataframe of reference db
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  seqnums_vector <- which(rn$n >= min_n)
  seq_range <- start_seq:length(seqnums_vector)
  
  loc_p <- c("Yeti yeti", "Sasquatch vulgaris")
  
  # define loop body once
  loop_body <- function(i) {
    if (verbose) message("Testing sequence ", i, " of ", length(seqnums_vector), ": ", rn[seqnums_vector[i], ]$s)
    tictoc::tic()
    
    # temporary patch to circumvent bug in blastg v0.1.4
    if (grepl("Orycteropus_afer", rn[seqnums_vector[i], ]$s)) {
      message("Skipping Aardvark...")
      tictoc::toc()
      return(invisible(NULL))
    }
    
    seq_num <- seqnums_vector[i]
    
    # diminished refdb
    refdb_dim <- LOSO_subset(refdb = r, seq_num, return_db = TRUE)
    refdb_looseq <- LOSO_subset(refdb = r, seq_num, return_db = FALSE)
    
    # unique tmp dir
    seq_tmp <- file.path(out, "tmp", paste0("seq_", i, "_", Sys.getpid()))
    dir.create(seq_tmp, recursive = TRUE, showWarnings = FALSE)
    
    # run ghostblaster
    ghost_data <- blastg(
      seqs = as.character(refdb_looseq),
      refdb = refdb_dim,
      out = seq_tmp,
      locals = loc_p,
      marker = marker,
      min_length = min_length,
      verbose = FALSE,
      BLAST_args = BLAST_args
    )
    
    write.csv(
      ghost_data$ghost_data,
      file = file.path(out, "loso", paste0("gd_", i, "_", ghost_data$ghost_data$qseqid[1], ".csv")),
      row.names = FALSE
    )
    
    # cleanup
    unlink(seq_tmp, recursive = TRUE)
    
    tictoc::toc()
    invisible(NULL)
  }
  
  # run either parallel or serial
  if (parallel) {
    stop("parallel=T disabled. Was skipping files.") # mclapply(seq_range, loop_body, mc.cores = n_cores)
  } else {
    for (i in seq_range) loop_body(i)
  }

}
