# function to extract a species from a reference database to allow leave-one-out benchmarking
  # inputs = reference database in dada2 format, species to leave out
  # output = diminished reference database OR sequences of only the "left out" species

LOOseq_refdb_subset <- function(refdb, # path to a reference fasta file (RDP format required)
                      seq_num, 
                      return_db = T,
                      verbose = F){

# load libraries
library(Biostrings, warn.conflicts = FALSE)
library(dplyr)
library(tidyr)

# inputs for testing
# refdb <- "data/midori_GB263_v16S.pga70.uni.l100.L300.n10.rdp.ama.fasta"
# seq_num <- 1 # sequence number ordered from 1 - total num sequences.
# return_db <- TRUE # if TRUE, return entire database minus LOO sequence; if FALSE, return just the LOO sequence
# verbose <- TRUE

# load reference database
if(class(refdb)[[1]] %in% "DNAStringSet"){
  r <- refdb}else{
    r <- readDNAStringSet(refdb)
  }

# get the fasta header information
header <- names(r)
  
# setup messages
if(verbose %in% TRUE){
message("Reading ", length(header), " sequences from reference database")
}

# extract single sequence
r1 <- r[seq_num]

# get header of extracted sequence
header.r1 <- names(r1)

# subtract LOO sequence from reference database
diminished_refdb <- r[grep(header.r1[1], header, invert = TRUE)]

# get species name and accession number of extracted sequence
r1.info <- data.frame(label = header.r1)  %>%
  separate(label, into = c("acc", "taxonomy"), sep = "\t") %>%
  separate(taxonomy,
           into = c("r", "k", "p", "c", "o", "f", "g", "s"),
           sep = ";") %>%
  mutate(species = gsub(" ", "_", s))

names(r1) <- paste0(r1.info$species, "_", r1.info$acc)

# print how many were extracted
if(verbose %in% TRUE){
message("Extracted ", length(header) - length(diminished_refdb),
        " sequence of ", r1.info$species)
if(return_db %in% TRUE){
  message("Returning: database without that sequence (n = ", length(diminished_refdb), ").")
}else{message("Returning: database of only that sequence.")}
}

# return desired output in DNAStringSet format
if(return_db %in% TRUE){
  return(diminished_refdb)
}else{return(r1)}

} # end of function
