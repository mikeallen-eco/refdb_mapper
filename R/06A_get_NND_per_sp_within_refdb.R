# get nearest-neighbor evolutionary (cophenetic) distance for each species within a reference database (NCBI)

library(Biostrings, warn.conflicts = FALSE)
library(dplyr)
library(tidyr)
library(ape)

get_NND_per_sp_within_refdb <- function(
             refdb_cur = refdb_cur_paths, # vector of paths to reference databases with NCBI species names
             min_n = 1, # only include species with ≥ min_n sequences
             tree = "data/phyl.tre",
             phyl_harmonized = phyl_harmonized_path,
             refdb_harmonized = refdb_harmonized_path,
             extinct = ncbi_extinct,
             verbose = TRUE) {
  
  refdb_harmonized <- read.csv(refdb_harmonized_path)
  phyl_harmonized <- read.csv(phyl_harmonized_path)
  
  nnd_list <- lapply(1:length(refdb_cur_paths), function(i){
    
  message("Processing refdb: \n", refdb_cur[i])

  # read in reference database
  if (grepl(refdb_cur[i], pattern = "fasta")) {
    r <- readDNAStringSet(refdb_cur_paths)
  } else{
    r <- refdb_cur[i]
  }
  
  # get df of reference database
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  
  # create taxonomy lookup
  species_list <- sort(unique(rn$s[which(rn$n >= min_n)])) %>%
    as.data.frame() %>%
    dplyr::rename(ncbi_name = 1) %>%
    filter(!grepl(ncbi_name, pattern = "_x_"),
           !ncbi_name %in% extinct) %>%
    left_join(refdb_harmonized %>% mutate(ncbi_name = gsub(" ", "_", full_sci_name)),
              by = join_by(ncbi_name)) %>%
    left_join(phyl_harmonized %>% 
                select(MOL_Accepted, phyl_name = full_sci_name) %>%
                group_by(MOL_Accepted) %>% slice_head(n = 1) %>% ungroup(),
              by = join_by(MOL_Accepted)) %>%
    mutate(mol_name = gsub(" ", "_", MOL_Accepted)) %>%
    select(ncbi_name, phyl_name, mol_name) %>%
    filter(!is.na(phyl_name))
  
  # reduce tree to species list
  reduced_tree <- drop.tip(phyl_tree, setdiff(phyl_tree$tip.label, underscore(unique(species_list$phyl_name))))
  
  # compute distance matrix
  dist.mat1 <- cophenetic.phylo(reduced_tree)
  row.names(dist.mat1) <- gsub("\\.", "-", row.names(dist.mat1)) # fix species with hyphens that got changed to periods
  colnames(dist.mat1) <- gsub("\\.", "-", colnames(dist.mat1)) # fix species with hyphens that got changed to periods
  
  # get order & sequence count from reference data
  rn_sp <- rn %>% 
    mutate(ncbi_name = ununderscore(s)) %>% 
    group_by(ncbi_name) %>% 
    summarize(n = length(ncbi_name),
              .groups = "drop")
  
 # get nearest neighbor co-phenetic distance for each phyl_name & add NCBI name
 NND_df <- apply(dist.mat1, 2, 
        FUN = function(x){
          sort(x)[2]
          }) %>%
   as.data.frame() %>%
   dplyr::rename(nnd = 1) %>%
   tibble::rownames_to_column("phyl_name") %>%
   mutate(phyl_name = ununderscore(phyl_name)) %>%
   left_join(species_list,
             by = join_by(phyl_name)) %>%
   left_join(rn_sp,
             by = join_by(ncbi_name)) %>%
   select(ncbi_name, phyl_name, nnd, n)
  
  return(NND_df)
  })
  
  return(nnd_list)
}
