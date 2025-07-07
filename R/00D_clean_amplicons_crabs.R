# function that takes a MIDORI2 RDP-format reference database and extracts amplicons

clean_amplicons_crabs <- function(input = "crabs_amplicons_pga.txt",
                                  out, l, L,
                                  db_name = "amplicon_source_date",
                                  taxon = "Mammalia",
                                  conda_dir = "/Users/mikea/miniconda3/bin/conda",
                                  conda_env = "crb2",
                                  verbose = FALSE) {

# dereplicate amplicon-species combinations
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--dereplicate", 
                 "--input", paste0(out, input),
                 "--output", paste0(out, "crabs_amplicons_pga.uni.txt"),
                 "--dereplication-method", "'unique_species'"), 
        stdout = TRUE, stderr = TRUE)

# quality filter 
# (length, number "N" bases, environmental samples, name issues like "sp.")
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--filter", 
                 "--input",
                 paste0(out, "crabs_amplicons_pga.uni.txt"),
                 "--output",
                 paste0(out, "crabs_amplicons_pga.uni.cln.txt"),
                 "--minimum-length", l,
                 "--maximum-length", L,
                 "--maximum-n", "10",
                 "--environmental --no-species-id --rank-na 3"), 
        stdout = TRUE, stderr = TRUE)

# export final RDP format fasta reference database 
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--export",
                 "--input",
                 paste0(out, "crabs_amplicons_pga.uni.cln.txt"), 
                 "--output", 
                 paste0(out, "refdb_", db_name, ".fasta"),
                 "--export-format 'rdp'"), 
        stdout = TRUE, stderr = TRUE)

# re-subset database by taxon
  # (required due to some higher taxonomy errors in midori2 that were fixed by importing into crabs)
r <- subset_raw_refdb_by_taxon(paste0(out, "refdb_", db_name, ".fasta"),
                               taxon = taxon)

# over-write final fasta with re-subsetted version
writeXStringSet(r, paste0(out, "refdb_", db_name, ".fasta"), append = FALSE)

# Remove temporary files created for crabs
unlink(paste0(out, "raw_seqs.fasta"))
unlink(paste0(out, "crabs_*"))

}