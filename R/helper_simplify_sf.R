library(rmapshaper)

simplify_sf <- function(original_rds = "~/Documents/mikedata/refdb_mapper/hydrobasin_map.rds", 
                        simple_rds = "data/hydrobasin_map_simple.rds",
                        p_keep = 0.05){
  
hydrobasins <- readRDS(original_rds)
hydrobasins_simple <- ms_simplify(hydrobasins, keep = 0.05, keep_shapes = TRUE)
saveRDS(hydrobasins_simple, simple_rds)
}