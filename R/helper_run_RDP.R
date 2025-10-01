# function to run the RDP taxonomic classifier in dada2 on a set of sequences 
# based on a RDP format reference database

# load libraries and helper functions
library(dada2)
library(tictoc)
library(dplyr)

run_RDP <- function(q_seqs, # named character vector of sequences (names = unique sequence IDs)
                    ref_seqs, # reference database as DNAStringSet in RDP format
                    out = "data/empirical/",
                    seed = 100){ # output directory
  
  tic()
  
  # make a temporary directory
  if (!dir.exists(paste0(out, "tmp/"))) {
    dir.create(paste0(out, "tmp/"))
  }
  
  # change header format to BayesANT
  r <- RDP_to_dada2(ref_seqs) # get dada2 format reference stringset
  rn <- RDP_to_dada2(ref_seqs, output = "df") # get dada2 format reference df for names
  
  # write formatted reference database as a fasta in a tmp folder
  writeXStringSet(r,
                  paste0(out,"tmp/tmp_RDP_refdb.fasta"),
                  append = F)
  
  set.seed(seed) # Initialize random number generator for reproducibility
  taxa.RDP <- assignTaxonomy(q_seqs, paste0(out, "tmp/tmp_RDP_refdb.fasta"), 
                             multithread=T, 
                             minBoot=0,
                             outputBootstraps = TRUE)
  
  # make a sequence<>qseqid lookup table
  ts_df <- data.frame(qseqid = names(q_seqs),
                      seqs = q_seqs)
  
  RDP.df <- taxa.RDP$tax %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    bind_cols(taxa.RDP$boot) %>% 
    dplyr::select(Kingdom = 2, Phylum = 3, Class = 4,
                  Order = 5, Family = 6, Genus = 7, Species = 8,
                  genboot = 14, spboot = 15, seqs = 1) %>%
    left_join(ts_df, by = join_by(seqs)) %>%
    suppressMessages()
  
  toc() 
  
  return(RDP.df)
  
}
