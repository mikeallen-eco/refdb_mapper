curate_amplicons <- function(refdb,
                             taxon = "Mammalia",
                             fwd,
                             rev,
                             out,
                             l,
                             L,
                             db_name,
                             dl_tax = FALSE,
                             conda_dir = "/Users/mikea/miniconda3/bin/conda",
                             conda_env = "crb2",
                             include_mismatched_primers = FALSE,
                             verbose = TRUE) {
  taxon_refdb <- subset_raw_refdb_by_taxon(refdb, taxon)
  
  if (dl_tax == TRUE) {
    download_NCBI_taxonomy_crabs(out = out)
  }
  
  extract_amplicons_crabs(
    ref_seqs = taxon_refdb,
    fwd = fwd,
    rev = rev,
    out = out,
    all_starts = include_mismatched_primers,
    conda_dir = conda_dir,
    conda_env = conda_env,
    verbose = verbose
  )
  
  clean_amplicons_crabs(input = "crabs_amplicons_pga.txt", out, l, L, db_name)
  
}