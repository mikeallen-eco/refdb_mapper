# function that takes a reference database and extracts amplicons using crabs
extract_amplicons_crabs <- function(ref_seqs = taxon_refdb, 
                                    fwd, rev, l, L,
                                    out,
                                    conda_dir = "/Users/mikea/miniconda3/bin/conda",
                                    conda_env = "crb2",
                                    crabs_path = "crabs", # or "python ~/path/to/reference_database_creator/crabs"
                                    tax_path = "/Users/mikea/Documents/mikedata/GhostBLASTer/BLASTrees/tax_20250609/",
                                    verbose = TRUE) {
  
  #--------------------------------------------------
  # Helper function to run system2 and save output
  #--------------------------------------------------
  run_and_log <- function(cmd, args, log_file) {
    cmd_out <- system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
    write(paste0("\n--- Command: ", paste(c(cmd, args), collapse = " "), "\n"), 
          file = log_file, append = TRUE)
    write(cmd_out, file = log_file, append = TRUE)
    invisible(cmd_out)
  }
  
  #--------------------------------------------------
  # Setup
  #--------------------------------------------------
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")  # e.g. "20251008_134205"
  log_file <- file.path(out, paste0("crabs_run_log_", timestamp, ".txt"))
  if (file.exists(log_file)) file.remove(log_file)  # start fresh
  
  imported_file <- file.path(out, "crabs_raw_seqs.txt")
  imported_already <- file.exists(imported_file)
  
  #--------------------------------------------------
  # Step 1: Import reference sequences into CRABS
  #--------------------------------------------------
  if (!imported_already) {
    if (verbose) message("Loading reference sequences into CRABS...")
    
    # detect if ref_seqs is a fasta file path or DNAStringSet
    if (grepl("fasta$", ref_seqs[1])) {
      ref_seqs_path <- ref_seqs
    } else {
      ref_seqs_path <- file.path(out, "raw_seqs.fasta")
      Biostrings::writeXStringSet(ref_seqs, ref_seqs_path)
    }
    
    if (verbose) tictoc::tic("Import complete...")
    
    import_args <- c(
      "run", "-n", conda_env, crabs_path,
      "--import",
      "--import-format", "NCBI",
      "--input", ref_seqs_path,
      "--names", file.path(tax_path, "names.dmp"),
      "--nodes", file.path(tax_path, "nodes.dmp"),
      "--acc2tax", file.path(tax_path, "nucl_gb.accession2taxid"),
      "--output", imported_file,
      "--ranks", "'superkingdom;phylum;class;order;family;genus;species'"
    )
    
    run_and_log(conda_dir, import_args, log_file)
    if (verbose) tictoc::toc()
    
  } else {
    message("crabs_raw_seqs.txt file detected. Skipping import step ...")
  }
  
  #--------------------------------------------------
  # Step 2: In silico PCR
  #--------------------------------------------------
  if (verbose) message("Performing in silico PCR ...")
  if (verbose) tictoc::tic("PCR complete ...")
  
  pcr_args <- c(
    "run", "-n", conda_env, crabs_path,
    "--in-silico-pcr",
    "--input", imported_file,
    "--output", file.path(out, "crabs_amplicons_pcr.txt"),
    "--forward", fwd,
    "--reverse", rev
  )
  
  run_and_log(conda_dir, pcr_args, log_file)
  if (verbose) tictoc::toc()
  
  #--------------------------------------------------
  # Step 3: Pairwise Global Alignment
  #--------------------------------------------------
  if (verbose) message("Performing PGA to extract amplicons without primer regions ...")
  if (verbose) tictoc::tic("PGA complete ...")
  
  pga_args <- c(
    "run", "-n", conda_env, crabs_path,
    "--pairwise-global-alignment",
    "--input", imported_file,
    "--amplicons", file.path(out, "crabs_amplicons_pcr.txt"),
    "--output", file.path(out, "crabs_amplicons_pga.txt"),
    "--forward", fwd,
    "--reverse", rev,
    "--size-select", "20000",
    "--percent-identity", "0.80",
    "--coverage", "95"
  )
  
  run_and_log(conda_dir, pga_args, log_file)
  if (verbose) tictoc::toc()
  
  #--------------------------------------------------
  # Done
  #--------------------------------------------------
  if (verbose) message("CRABS amplicon extraction complete. Log saved to: ", log_file)
  invisible(log_file)
}

