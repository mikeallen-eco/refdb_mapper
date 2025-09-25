library(Biostrings)
library(dplyr)
library(tictoc)
library(parallel)

LOSpO_ghostblaster <- function(refdb,
                               out,
                               marker = "average",
                               start_seq = 1,
                               blast_path = "/usr/local/ncbi/blast/bin",
                               BLAST_args = "-max_target_seqs 5000 -perc_identity 70",
                               min_length = 90,
                               verbose = TRUE,
                               parallel = FALSE,
                               n_cores = detectCores() - 1) {
  
  # create directories if needed
  if (!dir.exists(out)) dir.create(out)
  if (!dir.exists(file.path(out, "lospo"))) dir.create(file.path(out, "lospo"))
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
  sp_names <- sort(unique(rn$s))
  seq_range <- start_seq:length(sp_names)
  
  loc_p <- c("Yeti yeti", "Sasquatch vulgaris")
  
  # define loop body once
  loop_body <- function(i) {
    species_name <- sp_names[i]
    
    if (verbose) message("Testing species ", i, " of ", length(sp_names), ": ", species_name)
    tictoc::tic()
    
    # skip aardvark
    if (grepl("Orycteropus_afer", species_name)) {
      message("Skipping Aardvark...")
      tictoc::toc()
      return(invisible(NULL))
    }
    
    # diminished refdb (minus this species)
    refdb_dim <- LOSpO_subset(refdb = r, species = species_name, return_db = TRUE)
    # LOSpO target (just this species)
    refdb_loosp <- LOSpO_subset(refdb = r, species = species_name, return_db = FALSE)
    
    # unique tmp dir
    sp_tmp <- file.path(out, "tmp", paste0("sp_", i, "_", Sys.getpid()))
    dir.create(sp_tmp, recursive = TRUE, showWarnings = FALSE)
    
    # run ghostblaster
    ghost_data <- blastg(
      seqs = as.character(refdb_loosp),
      refdb = refdb_dim,
      out = sp_tmp,
      locals = loc_p,
      marker = marker,
      min_length = min_length,
      verbose = FALSE,
      BLAST_args = BLAST_args
    )
    
    # write results
    write.csv(
      ghost_data$ghost_data,
      file = file.path(out, "lospo", paste0("gd_lospo_", i, "_", species_name, ".csv")),
      row.names = FALSE
    )
    
    # cleanup
    unlink(sp_tmp, recursive = TRUE)
    
    tictoc::toc()
    invisible(NULL)
  }
  
  # run either parallel or serial
  if (parallel) {
    stop("paralkel=T disabled. Was skipping files.") # mclapply(seq_range, loop_body, mc.cores = n_cores)
  } else {
    for (i in seq_range) loop_body(i)
  }
}
