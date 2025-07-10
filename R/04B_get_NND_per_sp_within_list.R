# get nearest-neighbor evolutionary (cophenetic) distance for each species within a list of species (e.g., MOL)

# test_list = (read.csv("~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv") %>% 
#   filter(HYBAS_ID %in% "765062684"))$sciname

get_NND_per_sp_within_list <- function(
             sp_list, # character vector of species names
             tree = "data/phyl.tre",
             phyltax = "data/phyltax.csv",
             sp_list_tax = "data/geotax.csv",
             verbose = FALSE) {
  
  # load libraries and functions
  # library(Biostrings, warn.conflicts = FALSE)
  library(dplyr)
  library(tidyr)
  library(ape)
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  
  # read in phylogenetic tree taxonomy
  phyl_tax <- read.csv(phyltax) %>%
    mutate(seq_species = case_when(!is.na(ncbi_name) ~ ncbi_name,
                                   is.na(ncbi_name) & !is.na(gbif_name) ~ gbif_name,
                                   TRUE ~ phyl_name))
  
  # read in input species list taxonomy
  geo_tax <- read.csv(sp_list_tax) %>%
    dplyr::rename(geo_name = orig_name) %>%
    mutate(seq_species = case_when(!is.na(ncbi_name) ~ ncbi_name,
                                   is.na(ncbi_name) & !is.na(gbif_name) ~ gbif_name,
                                   TRUE ~ geo_name))
  
  # get list of species (NCBI) with ≥ min_n variants present that can be matched to phylogeny
    # note: 683 of 724 mammals were matchable, the remainder mainly extinct or domestic species
  species_list <- ununderscore(unique(sp_list)) %>%
    as.data.frame() %>%
    dplyr::rename(geo_name = 1) %>%
    left_join(geo_tax, by = join_by(geo_name)) %>%
    left_join(phyl_tax %>% select(seq_species, phyl_name), 
              by = join_by(seq_species)) %>%
    group_by(seq_species) %>%
    slice_head(n = 1) %>%
    mutate(phyl_name = case_when(is.na(phyl_name) ~ seq_species,
                                 TRUE ~ phyl_name),
           in_phyl = case_when(seq_species %in% phyl_tax$seq_species ~ 1,
                               TRUE ~ 0)) %>%
    select(order, geo_name, seq_species, ncbi_name, gbif_name, phyl_name, in_phyl) %>%
    filter(!is.na(phyl_name)) %>%
    group_by(seq_species) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # reduce tree to species list
  reduced_tree_names <- unique(gsub(" ", "_", species_list$phyl_name))
  reduced_tree <- drop.tip(phyl_tree, setdiff(phyl_tree$tip.label, reduced_tree_names))
  
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
   full_join(species_list,
             by = join_by(phyl_name)) %>%
   select(order, geo_name, seq_species, ncbi_name, phyl_name, in_phyl, nnd) %>%
   arrange(order, geo_name)
  
  return(NND_df)
  
}