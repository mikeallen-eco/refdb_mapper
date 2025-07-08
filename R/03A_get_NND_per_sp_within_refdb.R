# get nearest-neighbor evolutionary (cophenetic) distance for each species within a reference database (NCBI)

get_NND_per_sp_within_refdb <- function(
             refdb, # reference database with NCBI species names
             min_n, # only include species with ≥ min_n sequences
             tree = "data/phyl.tre",
             phyltax = "data/phyltax.csv",
             verbose = FALSE) {
  
  # load libraries and functions
  library(Biostrings, warn.conflicts = FALSE)
  library(dplyr)
  library(tidyr)
  library(ape)
  
  # read in reference database
  if (grepl(refdb[1], pattern = "fasta")) {
    r <- readDNAStringSet(refdb)
  } else{
    r <- refdb
  }
  
  # get df of reference database
  rn <- RDP_to_dataframe(r) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup()
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  
  # read in phylogenetic tree taxonomy
  phyl_tax <- read.csv(phyltax) %>%
    mutate(seq_species = case_when(!is.na(ncbi_name) ~ ncbi_name,
                                   is.na(ncbi_name) & !is.na(gbif_name) ~ gbif_name,
                                   TRUE ~ phyl_name))
  
  # get list of species (NCBI) with ≥ min_n variants present that can be matched to phylogeny
    # note: 683 of 724 mammals were matchable, the remainder mainly extinct or domestic species
  species_list <- gsub("_", " ", sort(unique(rn$s[which(rn$n > min_n)]))) %>%
    as.data.frame() %>%
    dplyr::rename(ncbi_name = 1) %>%
    left_join(phyl_tax %>% select(ncbi_name, phyl_name_via_ncbi = phyl_name), by = join_by(ncbi_name)) %>%
    left_join(phyl_tax %>% select(ncbi_name = gbif_name,
                                  phyl_name_via_gbif = phyl_name) %>%
                group_by(ncbi_name) %>% slice_head(n = 1), by = join_by(ncbi_name)) %>%
    left_join(phyl_tax %>% select(ncbi_name = phyl_name) %>% mutate(phyl_name_via_phyl = ncbi_name) %>%
                group_by(ncbi_name) %>% slice_head(n = 1), by = join_by(ncbi_name)) %>%
    mutate(phyl_name = case_when(!is.na(phyl_name_via_ncbi) ~ phyl_name_via_ncbi,
                                 is.na(phyl_name_via_ncbi) & 
                                   !is.na(phyl_name_via_gbif) ~ phyl_name_via_gbif,
                                 TRUE ~ phyl_name_via_phyl)) %>%
    select(ncbi_name, phyl_name) %>%
    filter(!is.na(phyl_name)) %>%
    group_by(ncbi_name) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # reduce tree to species list
  reduced_tree <- drop.tip(phyl_tree, setdiff(phyl_tree$tip.label, underscore(species_list$phyl_name)))
  
  # compute distance matrix
  dist.mat1 <- cophenetic.phylo(reduced_tree)
  row.names(dist.mat1) <- gsub("\\.", "-", row.names(dist.mat1)) # fix species with hyphens that got changed to periods
  colnames(dist.mat1) <- gsub("\\.", "-", colnames(dist.mat1)) # fix species with hyphens that got changed to periods
  
 # get nearest neighbor co-phenetic distance for each phyl_name & add NCBI name
 NND_df <- apply(dist.mat1, 2, 
        FUN = function(x){
          sort(x)[2]
          }) %>%
   as.data.frame() %>%
   dplyr::rename(nnd = 1) %>%
   tibble::rownames_to_column("phyl_name") %>%
   mutate(phyl_name = ununderscore(phyl_name)) %>%
   left_join(species_list_ncbi,
             by = join_by(phyl_name)) %>%
   left_join(rn %>% mutate(ncbi_name = ununderscore(s)) %>% select(ncbi_name, order = o) %>% distinct(),
             by = join_by(ncbi_name)) %>%
   select(ncbi_name, phyl_name, order, nnd)
  
  return(NND_df)
  
}