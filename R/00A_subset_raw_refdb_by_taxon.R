# subset a reference database in fasta format by a character string in the header (e.g., taxon name)

library(Biostrings, warn.conflicts = FALSE)

subset_raw_refdb_by_taxon <- function(refdb,
                                      taxon){ 

refdb_seqs <- readDNAStringSet(refdb)
headers <- names(refdb_seqs)
refdb_seqs_taxon <- refdb_seqs[grep(taxon, headers, invert = FALSE)]

return(refdb_seqs_taxon)

}