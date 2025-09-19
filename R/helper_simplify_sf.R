library(rmapshaper)

simplify_sf <- function(original_map = hydrobasin_map, 
                        p_keep = 0.05){
  
hydrobasins_simple <- ms_simplify(original_map, keep = 0.05, keep_shapes = TRUE)
hydrobasins_simple
}