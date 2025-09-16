curate_amplicons <- function(refdb = raw_refdb_path,
                             taxon = "Mammalia",
                             fwd = fwd,
                             rev = rev,
                             out = out_path,
                             l = l,
                             L = L,
                             db_name = db_name,
                             dl_tax = FALSE,
                             conda_dir = "/Users/mikea/miniconda3/bin/conda",
                             conda_env = "crb2",
                             crabs_path = "crabs", # or e.g., for newest version "python ~/reference_database_creator/crabs"
                             keep_imported = TRUE,
                             verbose = TRUE) {
  
  taxon_refdb <- subset_raw_refdb_by_taxon(refdb, taxon)
  
  if (dl_tax == TRUE) {
    download_NCBI_taxonomy_crabs(out = out,
                                 conda_dir = conda_dir,
                                 conda_env = conda_env,
                                 crabs_path = crabs_path)
  }
  
  extract_amplicons_crabs(
    ref_seqs = taxon_refdb,
    fwd = fwd,
    rev = rev,
    out = out,
    conda_dir = conda_dir,
    conda_env = conda_env,
    crabs_path = crabs_path,
    verbose = verbose
  )
  
  clean_amplicons_crabs(input = "crabs_amplicons_pga.txt", out, l, L, db_name,
                        conda_dir = conda_dir,
                        conda_env = conda_env,
                        crabs_path = crabs_path)
  
}