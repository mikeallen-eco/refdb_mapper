# get nearest-neighbor evolutionary (cophenetic) distance for each species within a reference database (NCBI)

library(Biostrings, warn.conflicts = FALSE)
library(dplyr)
library(tidyr)
library(ape)

get_NND_per_sp_within_refdb <- function(
             refdb_cur = refdb_cur_paths, # vector of paths to reference databases with NCBI species names
             tree = "data/phyl.tre",
             refdb_harmonized_path = refdb_harmonized_path,
             mol_to_phyl_harmonized = mol_to_phyl_harmonized_path,
             extinct = c(ncbi_extinct, phyl_extinct)) {
  
  nnd_list <- lapply(1:length(refdb_cur_paths), function(i){
    
  message("Processing refdb: \n", refdb_cur[i])

  # read in reference database
  if (grepl(refdb_cur[i], pattern = "fasta")) {
    r <- readDNAStringSet(refdb_cur_paths[[i]])
  } else{
    r <- refdb_cur[i]
  }
  
    refdb_harmonized <- read.csv(refdb_harmonized_path)

  # get df of reference database
  rn <- RDP_to_dataframe(r) %>%
    select(ncbi_name = s) %>%
    distinct() %>%
    left_join(refdb_harmonized %>% mutate(ncbi_name = underscore(full_sci_name)),
              by = join_by(ncbi_name)) %>%
    filter(!grepl(ncbi_name, pattern = "_x_"),
           !ncbi_name %in% extinct) %>%
    select(-full_sci_name) %>%
    filter(!is.na(BB_Accepted))
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  
 # get nearest neighbor co-phenetic distance for each phyl_name & add NCBI name
 NND_df <- get_NND_per_sp_within_list(sp_list = rn$BB_Accepted, 
                                      tree = tree, 
                                      mol_to_phyl_harmonized_path = mol_to_phyl_harmonized,
                                      extinct = extinct)
  
  return(NND_df)
  })
  
  return(nnd_list)
}
