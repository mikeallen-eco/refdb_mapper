# function to reformat a DNAStringSet from RDP format to dada2 taxonomy assignment format

RDP_to_dada2 <- function(rdp_DNAStringSet,
                         output = "stringset"){ # or "df"
  
  # format into dada2 fasta header
  parsed_df <- data.frame(header = names(rdp_DNAStringSet)) %>%
    separate(header, into = c("id", "taxonomy"), sep = "\t") %>%
    separate(taxonomy, into = c("r", "k", "p", "c", "o", "f", "g", "s"), 
             sep = ";") %>%
    mutate(hdr = paste(k,p,c,o,f,g,s, sep = ";")) 
  
  # check_taxonomy_consistency(df = parsed_df)
  
  # add BayesANT header to original DNA String Set
  reformatted_DNAStringSet <- rdp_DNAStringSet
  names(reformatted_DNAStringSet) <- parsed_df$hdr
  
  if(output %in% "stringset"){
    return(reformatted_DNAStringSet)}else{
      return(parsed_df)
    }
  
}