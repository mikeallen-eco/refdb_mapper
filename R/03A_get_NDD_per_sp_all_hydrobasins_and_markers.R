# get evolutionary nearest neighbor distance for species within all hydrobasins for each marker

# load libraries and functions
library(dplyr)
library(tidyr)
library(ape)
library(purrr)

get_NDD_per_sp_all_hydrobasins_and_markers <- function(hydrobasin_refdb_info = hydrobasin_refdb_info,
                                                       markers = markers,
                                                       n_cores = 4,
                                                       tree_path = "data/phyl.tre",
                                                       mol_to_phyl_path = mol_to_phyl_harmonized_path,
                                                       extinct = c(ncbi_extinct, phyl_extinct)) {
  # read in the phylogenetic tree (with temp patch to insert humans)
  phyl_tree <- ape::read.tree(tree_path) # combined tree covering all vertebrates
  phyl_tree$tip.label[grepl(phyl_tree$tip.label, pattern = "Homo_neanderthalensis")] <- "Homo_sapiens"
  
  # set up MOL -> phylogeny lookup information
  mol_phyl <- read.csv(mol_to_phyl_path) %>%
    mutate(mol_name = underscore(full_sci_name),
           phyl_name = underscore(BB_Accepted)) %>%
    filter(!phyl_name %in% extinct, !is.na(phyl_name))
  
  # reduce tree to species list
  reduced_tree <- drop.tip(phyl_tree, setdiff(phyl_tree$tip.label, unique(mol_phyl$phyl_name)))
  
  # compute distance matrix
  dist_full <- cophenetic.phylo(reduced_tree)
  row.names(dist_full) <- gsub("\\.", "-", row.names(dist_full)) # fix species with hyphens that got changed to periods
  colnames(dist_full) <- gsub("\\.", "-", colnames(dist_full)) # fix species with hyphens that got changed to periods
  
  hybas_ids <- unique(hydrobasin_refdb_info$HYBAS_ID)
  
  hybas_nnd_by_marker <- list()
  
  for (i in seq_along(markers)) {
    message("Performing NND calculations by hydrobasin for marker: ",
            markers[i])
    nnd_df <- mclapply(seq_along(hybas_ids), function(y) {
      message("Processing hybas_id ", y, " of ", length(hybas_ids))
      get_NND_for_hydrobasin(
        hydrobasin_refdb_info,
        dist_full,
        mol_phyl,
        marker = markers[i],
        HYBAS_ID_num = hybas_ids[y]
      )
      }, mc.cores = n_cores) %>%
      do.call(bind_rows, .)
    
    hybas_nnd_by_marker[[i]] <- nnd_df
    
  }
  
  # Join all marker-specific dfs together
  hybas_nnd_all <- purrr::reduce(hybas_nnd_by_marker,
                                 full_join,
                                 by = c("HYBAS_ID", "mol_name", "phyl_name")) %>%
    left_join(hydrobasin_refdb_info %>% select(HYBAS_ID, mol_name, order, family),
              by = join_by(mol_name, HYBAS_ID)) %>%
    dplyr::select(HYBAS_ID,
                  order,
                  family,
                  mol_name,
                  phyl_name,
                  contains("12S"),
                  contains("16S"))
  
  return(hybas_nnd_all)
  
}
