# get nearest-neighbor evolutionary (cophenetic) distance for each species within a list of species (e.g., MOL)

# test_list = (read.csv("~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv") %>% 
#   filter(HYBAS_ID %in% "765062684"))$sciname

get_NND_per_sp_within_list <- function(
             sp_list, # character vector of species names
             tree = "data/phyl.tre",
             phyltax = "data/phyltax.csv",
             sp_list_tax = "data/geotax.csv",
             verbose = FALSE) {
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  


  if(length(sp_list) > 1){ 
  # reduce tree to species list
  reduced_tree_names <- unique(gsub(" ", "_", unique(species_list$phyl_name)))
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
   select(order, geo_name, ncbi_name, phyl_name, in_phyl, nnd) %>%
   arrange(order, geo_name)
  }else{NND_df <- species_list %>%
    select(order, geo_name, ncbi_name, phyl_name, in_phyl) %>%
    mutate(nnd = NA)}
  
  return(NND_df)
  
}