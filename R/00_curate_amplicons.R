curate_amplicons <- function(refdb,
                             taxon,
                             fwd,
                             rev,
                             out,
                             l,
                             L,
                             db_name,
                             conda_dir,
                             conda_env,
                             verbose = FALSE) {
  taxon_refdb <- subset_raw_refdb_by_taxon(refdb, taxon)
  
  if (dl_tax == TRUE) {
    download_NCBI_taxonomy_crabs(out = out)
  }
  
  extract_amplicons_crabs(
    ref_seqs = taxon_refdb,
    fwd = fwd,
    rev = rev,
    out = out,
    verbose = TRUE
  )
  
  clean_amplicons_crabs(input = "crabs_amplicons_pga.txt", out, l, L, db_name)
  
}