# get evolutionary nearest neighbor distance of each species within all hydrobasins in a gpkg

library(data.table)

get_NDD_per_sp_all_hydrobasins <- function(hydrobasin_map,
                                           hydrobasin_species = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv",
                                           tree_path = "data/phyl.tre",
                                           phyltax_path = "data/phyltax.csv",
                                           sp_list_tax_path = "data/geotax.csv",
                                           verbose = T){
  
  # read in map of hydrobasins
  if (class(hydrobasin_map)[1] %in% "character") {
    hydrobasins <- st_read(hydrobasin_map)
  } else{
    hydrobasins <- hydrobasin_map
  }
  
  hydrosp_all <- fread(hydrobasin_species)
  
  hybas_ndd_list <- lapply(1:length(hydrobasins$HYBAS_ID), function(x){
    print(x)
    
    hybas_sp_list <- hydrosp_all %>%
      filter(HYBAS_ID %in% hydrobasins$HYBAS_ID[x])
    
    df <- get_NND_per_sp_within_list(sp_list = hybas_sp_list$sciname,
                               tree = tree_path,
                               phyltax = phyltax_path,
                               sp_list_tax = sp_list_tax_path,
                               verbose = verbose) %>%
      mutate(HYBAS_ID = hydrobasins$HYBAS_ID[x])
    return(df)
  })
  
  return(hybas_ndd_list)
  
}
