# subset a reference database in fasta format by a character string in the header (e.g., taxon name)

library(Biostrings, warn.conflicts = FALSE)
library(yaml, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(argparse)

subset_raw_refdb_by_taxon <- function(refdb,
                                      taxon = "Mammalia",
                                      output_file){ 

refdb_seqs <- readDNAStringSet(refdb)
headers <- names(refdb_seqs)
refdb_seqs_taxon <- refdb_seqs[grep(taxon, headers, invert = FALSE)]

Biostrings::writeXStringSet(refdb_seqs_taxon, output_file)

}

# Inputs from command line
# input_file, primer_yaml, output_path_linked, output_path_unlinked, cutadapt_error_rate, n_cores

# Define the argument parser
parser <- ArgumentParser()

# Add named arguments
parser$add_argument("--refdb", type = "character", help = "input fasta file to be subset")
parser$add_argument("--taxon", type = "character", help = "Taxon name to be extracted - e.g., Mammalia")
parser$add_argument("--output_file", type = "character", help = "output path for subset fasta")

# Parse the command line arguments
args <- parser$parse_args()

# actual call to function, run_cutadapt
subset_raw_refdb_by_taxon(refdb = args$refdb, 
                          taxon = args$taxon,
                          output_file = args$output_file)
