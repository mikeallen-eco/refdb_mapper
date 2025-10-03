# get global watershed attributes for a set of coordinates

library(sf)
library(dplyr)

 get_polygon_attributes_from_coords <- function(lon = -59.912786,
                                                lat = -2.069209,
                                                include_marine = FALSE,
                                                polygons = final_sf,
                                                verbose = FALSE) {
  # temporarily turn off the S2 engine
  sf::sf_use_s2(FALSE) %>% suppressMessages()
   
  # define marine families
   marine <- c("Balaenopteridae", "Delphinidae", "Ziphiidae", "Dugongidae", 
               "Kogiidae", "Physeteridae", "Balaenidae", "Cetotheriidae", 
               "Otariidae", "Trichechidae", "Phocidae", "Phocoenidae", 
               "Monodontidae", "Odobenidae")
  
  # Make point sf
  pt <- st_sf(
    geometry = st_sfc(st_point(c(lon, lat)), crs = 4326)
  )
  
  # filter out marine species if applicable
  if (include_marine == FALSE) { polygons <- polygons %>% filter(!family %in% marine)}
  
  # Ensure both are in the same CRS
  if (st_crs(polygons) != st_crs(pt)) {
    pt <- st_transform(pt, st_crs(polygons))
  }

  # Spatial join: find the polygon that contains the point
  if(verbose) message("Getting attributes from polygon data ...")
  result <- st_join(pt, polygons, join = st_intersects) %>% suppressMessages()
  if(verbose) message("Getting attributes from polygon data ...")
  
  # Drop geometry to just return attributes
  result_df <- result %>% st_drop_geometry()
  
  # turn S2 engine back on when done
  sf::sf_use_s2(TRUE) %>% suppressMessages()
  
  return(result_df)
}
