# get nearest-neighbor evolutionary (cophenetic) distance for each species within a reference database (NCBI)

library(Biostrings, warn.conflicts = FALSE)
library(dplyr)
library(tidyr)
library(ape)

build_error_model_data <- function(
             refdb_cur = refdb_cur_paths, # vector of paths to reference databases with NCBI species names
             tree = "data/phyl.tre",
             refdb_harmonized_path = refdb_harmonized_path,
             mol_to_phyl_harmonized = mol_to_phyl_harmonized_path,
             extinct = c(ncbi_extinct, phyl_extinct),
             marker_names = markers) {
  
  nnd_list <- lapply(1:length(refdb_cur_paths), function(i){
    
  message("Building predictor data for marker: ", marker_names[i])

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
    left_join(refdb_harmonized %>% mutate(ncbi_name = underscore(full_sci_name)),
              by = join_by(ncbi_name)) %>%
    filter(!grepl(ncbi_name, pattern = "_x_"),
           !ncbi_name %in% extinct) %>%
    select(mol_name = BB_Accepted) %>%
    filter(!is.na(mol_name)) %>%
    group_by(mol_name) %>%
    summarise(n_seqs = length(mol_name)) %>%
    mutate(mol_name = underscore(mol_name),
           marker = marker_names[i])
  
 # get nearest neighbor co-phenetic distance for each phyl_name & add NCBI name
  NND_df <- get_NND_per_sp_within_list(
    sp_list = rn$mol_name,
    tree = tree,
    mol_to_phyl_harmonized_path = mol_to_phyl_harmonized,
    extinct = extinct
  ) %>%
    left_join(rn, by = join_by(mol_name)) %>%
    mutate(across(.cols = c(contains("_12S"), contains("_16S")), 
                  ~ replace_na(.x, 0)))
  
  return(NND_df)
  })
  
  nnd_df <- do.call(bind_rows, nnd_list)
  return(nnd_df)
}
