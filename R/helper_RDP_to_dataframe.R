# function to convert a RDP format DNAStringSet into a data frame

library(dplyr)
library(Biostrings)
library(tidyr)

RDP_to_dataframe <- function(RDPstringset,
                             include_seqs = F){
  
  df <- data.frame(header = names(RDPstringset)) %>%
    separate(header, into = c("id", "taxonomy"), sep = "\t") %>%
    separate(taxonomy, into = c("r", "k", "p", "c", "o", "f", "g", "s"), 
             sep = ";")
  
  if(include_seqs %in% TRUE) {
    df$seqs <- as.character(RDPstringset)
  }
  
  return(df)

  }