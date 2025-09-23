# get evolutionary nearest neighbor distance of each species within all hydrobasins in a gpkg

# load libraries and functions
library(data.table)
library(dplyr)
library(tidyr)
library(ape)

get_NDD_per_sp_all_hydrobasins <- function(hydrobasin_map = hydrobasin_map,
                                           hydrobasin_species = hydrobasin_species_path,
                                           tree_path = "data/phyl.tre",
                                           mol_to_phyl_harmonized = mol_to_phyl_harmonized_path,
                                           extinct = c(ncbi_extinct, phyl_extinct)){
  
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }
  
  hydrosp_all <- fread(hydrobasin_species) %>%
    mutate(mol_name = underscore(sciname))
  n_hybas <- length(hydrobasins$HYBAS_ID)
  
  hybas_nnd_list <- lapply(1:length(hydrobasins$HYBAS_ID), function(x){
    message(x, " of ", n_hybas)
    
    hybas_sp_list <- hydrosp_all %>%
      filter(HYBAS_ID %in% hydrobasins$HYBAS_ID[x])
    
    if(nrow(hybas_sp_list) > 0){
    df <- get_NND_per_sp_within_list(sp_list = hybas_sp_list$mol_name,
                               tree = tree_path,
                               mol_to_phyl_harmonized_path = mol_to_phyl_harmonized,
                               extinct = extinct) %>%
      right_join(hybas_sp_list %>% select(mol_name), by = join_by(mol_name)) %>%
      mutate(HYBAS_ID = hydrobasins$HYBAS_ID[x])
    }else{df <- data.frame(HYBAS_ID = hydrobasins$HYBAS_ID[x])}
    
    return(df)
  })
  
  return(hybas_nnd_list)
  
}
