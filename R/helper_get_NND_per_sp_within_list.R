# get nearest-neighbor evolutionary (cophenetic) distance for each species within a list of species (e.g., MOL)

# test_list = (read.csv("~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv") %>% 
#   filter(HYBAS_ID %in% "765062684"))$sciname

get_NND_per_sp_within_list <- function(
             sp_list = unique(rn$BB_Accepted), # character vector of MOL species names
             tree = "data/phyl.tre",
             mol_to_phyl_harmonized_path = mol_to_phyl_harmonized_path,
             extinct = c(ncbi_extinct, phyl_extinct)) {
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates
  
  # format MOL species list
  sp_list <- sort(unique(underscore(sp_list)))
  
  # create taxonomy lookup for MOL -> phylogeny name
  mol_phyl <- read.csv(mol_to_phyl_harmonized_path) %>%
    select(mol_name = uid, phyl_name = BB_Accepted) %>%
    mutate(phyl_name = underscore(phyl_name)) %>%
    filter(mol_name %in% sp_list, 
           !phyl_name %in% extinct) %>%
    filter(!is.na(phyl_name))
  
  # reduce tree to species list
  if(length(mol_phyl$mol_name) > 1){ 
  reduced_tree <- drop.tip(phyl_tree, setdiff(phyl_tree$tip.label, mol_phyl$phyl_name))
  
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
    left_join(mol_phyl,
              by = join_by(phyl_name)) %>%
    select(mol_name, phyl_name, nnd) %>%
    arrange(mol_name)
  }else{
      NND_df <- mol_phyl %>%
  select(mol_name, phyl_name) %>%
  mutate(nnd = NA) %>%
        arrange(mol_name)}
  
  dup_phyl <- unique(NND_df$phyl_name[duplicated(NND_df$phyl_name)])
  
  NND_df_final <- NND_df %>%
    mutate(nnd = case_when(phyl_name %in% dup_phyl ~ 0,
                           TRUE ~ nnd))
  
  return(NND_df_final)
  
}