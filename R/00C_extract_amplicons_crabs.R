# function that takes a MIDORI2 RDP-format reference database and extracts amplicons
extract_amplicons_crabs <- function(ref_seqs = taxon_refdb, 
                                    fwd, rev, l, L,
                                    out,
                                    conda_dir = "/Users/mikea/miniconda3/bin/conda",
                                    conda_env = "crb2",
                                    crabs_path = "crabs", # or "python ~/path/to/reference_database_creator/crabs"
                                    tax_path = "/Users/mikea/Documents/mikedata/GhostBLASTer/BLASTrees/tax_20250609/",
                                    verbose = TRUE) {
  
  #--------------------------------------------------
  # Helper function to run system2 and log output
  #--------------------------------------------------
  run_and_log <- function(cmd, args, log_file) {
    log_con <- file(log_file, open = "a")
    on.exit(close(log_con))
    system2(cmd, args = args, stdout = log_con, stderr = log_con)
  }
  
  #--------------------------------------------------
  # Setup
  #--------------------------------------------------
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
  }
  
  log_file <- file.path(out, "crabs_run.log")
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
    "--coverage", "0.95"
  )
  
  run_and_log(conda_dir, pga_args, log_file)
  if (verbose) tictoc::toc()
  
  #--------------------------------------------------
  # Done
  #--------------------------------------------------
  if (verbose) message("CRABS amplicon extraction complete. Log saved to: ", log_file)
  invisible(log_file)
}


# # function that takes a MIDORI2 RDP-format reference database and extracts amplicons
# 
# extract_amplicons_crabs <- function(ref_seqs = taxon_refdb, 
#                                     fwd, rev, l, L,
#                                    out,
#                                    conda_dir = "/Users/mikea/miniconda3/bin/conda",
#                                    conda_env = "crb2",
#                                    crabs_path = "crabs", # or "python ~/path/to/reference_database_creator/crabs" to use the newest version
#                                    tax_path = "/Users/mikea/Documents/mikedata/GhostBLASTer/BLASTrees/tax_20250609/",
#                                    verbose = TRUE){
#   
#   # make working directory if needed
#   if (!dir.exists(out)) {
#     dir.create(out)
#   }
#   
#   # detect if imported crabs format exists already
#   imported_already <- file.exists(paste0(out, "crabs_raw_seqs.txt"))
#   
#   if(imported_already == FALSE){
#     
#   if (verbose == TRUE) {message("Loading reference sequences into CRABS...")}
#   # detect if refseqs is a fasta file and if not treat it like a DNAStringset
#   if(grepl(ref_seqs[1], pattern = "fasta")){
#     ref_seqs_path <- ref_seqs
#   }else{
#     writeXStringSet(ref_seqs, paste0(out, "raw_seqs.fasta"))
#     ref_seqs_path <- paste0(out, "raw_seqs.fasta")
#   }
# 
#   if (verbose == TRUE) {tictoc::tic("Import complete...")}
# system2(conda_dir, 
#         args = c("run", "-n", conda_env, crabs_path, 
#                  "--import",
#                  "--import-format", "NCBI",
#                  "--input", ref_seqs_path,
#                  "--names", paste0(tax_path, "names.dmp"),
#                  "--nodes", paste0(tax_path, "nodes.dmp"),
#                  "--acc2tax", paste0(tax_path, "nucl_gb.accession2taxid"),
#                  "--output", paste0(out, "crabs_raw_seqs.txt"),
#                  "--ranks", "'superkingdom;phylum;class;order;family;genus;species'"
#         ), 
#         stdout = TRUE, stderr = TRUE)
# if (verbose == TRUE) {tictoc::toc()}
# }else{message("crabs_raw_seqs.txt file detected. Skipping import step ...")}
# 
# if (verbose == TRUE) {message("Performing in silico PCR ...")}
# # In silico PCR to extract amplicons with primer regions (w/ ≤ 4 mismatches)
# if (verbose == TRUE) {tictoc::tic("PCR complete ...")}
# system2(conda_dir, 
#         args = c("run", "-n", conda_env, crabs_path, 
#                  "--in-silico-pcr",
#                  "--input", paste0(out, "crabs_raw_seqs.txt"),
#                  "--output", paste0(out, "crabs_amplicons_pcr.txt"),
#                  "--forward", fwd,
#                  "--reverse", rev), 
#         stdout = TRUE, stderr = TRUE)
# if (verbose == TRUE) {tictoc::toc()}
# 
# if (verbose == TRUE) {message("Performing pga to extract amplicons without primer regions ...")}
# if (verbose == TRUE) {tictoc::tic("pga complete ...")}
# pga_arg_vector <- c("run", "-n", conda_env, crabs_path, 
#                     "--pairwise-global-alignment",
#                     "--input", paste0(out, "crabs_raw_seqs.txt"),
#                     "--amplicons", paste0(out, "crabs_amplicons_pcr.txt"),
#                     "--output", paste0(out, "crabs_amplicons_pga.txt"),
#                     "--forward", fwd,
#                     "--reverse", rev,
#                     "--size-select", "20000", 
#                     "--percent-identity", "0.70",
#                     "--coverage", "0.95")
# system2(conda_dir, 
#         args = pga_arg_vector, 
#         stdout = TRUE, stderr = TRUE)
# if (verbose == TRUE) {tictoc::toc()}
# 
# }
