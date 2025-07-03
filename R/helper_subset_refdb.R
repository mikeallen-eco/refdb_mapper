# function to subset a reference database by a MOL species list
  # requires reference database in DNAStringSet format with RDP header
  # locals is a character vector of MOL names

library(dplyr)
library(Biostrings)
library(tidyr)

subset_refdb <- function(RDPstringset,
                         locals){
  
  # load in taxonomy look up table
  loc_tax <- read.csv("data/geotax.csv") %>%
    filter(orig_name %in% locals) %>%
    mutate(orig_name = gsub(" ", "_", orig_name),
           ncbi_name = gsub(" ", "_", ncbi_name),
           gbif_name = gsub(" ", "_", gbif_name))  
  
  # read in string set, parse header, and subset by MOL names list
  loc_refdb_df <- data.frame(header = names(RDPstringset)) %>%
    separate(header, into = c("id", "taxonomy"), sep = "\t") %>%
    separate(taxonomy, into = c("r", "k", "p", "c", "o", "f", "g", "s"), 
             sep = ";") %>%
    mutate(seqs = as.character(RDPstringset),
           hdr = paste0(id, "\t", r, ";", k, ";", p, ";", c, ";", 
                        o, ";", f, ";", g, ";", s)) %>%
    filter(s %in% loc_tax$orig_name |
             s %in% loc_tax$ncbi_name |
             s %in% loc_tax$gbif_name)
  
  # create subset DNAStringSet
  loc_refdb <- DNAStringSet(loc_refdb_df$seqs)
  names(loc_refdb) <- loc_refdb_df$hdr
  
  return(loc_refdb)
  
}