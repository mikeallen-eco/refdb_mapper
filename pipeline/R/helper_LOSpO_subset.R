# function to extract a species from a reference database to allow leave-one-out benchmarking
  # inputs = reference database in dada2 format, species to leave out
  # output = diminished reference database OR sequences of only the "left out" species

LOSpO_subset <- function(refdb, 
                      species, 
                      return_db = T,
                      verbose = F){

# load libraries
library(Biostrings, warn.conflicts = FALSE)
library(dplyr)
library(tidyr)

# load reference database
if(grepl(refdb, pattern = "fasta")[1]){
    r <- readDNAStringSet(refdb)
  }else{r <- refdb}

# get the fasta header information
header <- names(r)

# format species name
sp <- gsub(" ", "_", species)
  
# setup messages
if(verbose %in% TRUE){
message("Reading reference database: ", 
        length(header), " sequences from ", 
        length(unique(header)), " species.")
}

# subtract LOO species from reference database
diminished_refdb <- r[grep(sp, header, invert = TRUE)]

r1 <- r[grep(sp, header, invert = FALSE)]
header.r1 <- names(r1)

# get species name and accession number of test sequences to add to header
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
        " sequences of ", species)
if(return_db %in% TRUE){
  message("Returning: database without those sequences.")
}else{message("Returning: database of only those sequences.")}
}

# return desired output in DNAStringSet format
if(return_db %in% TRUE){
  return(diminished_refdb)
}else{return(r1)}

} # end of function
