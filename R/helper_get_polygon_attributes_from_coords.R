# get global watershed attributes for a set of coordinates

library(sf)
library(dplyr)

 get_polygon_attributes_from_coords <- function(lon = -59.912786,
                                                lat = -2.069209,
                                                polygons = final_sf) {
  # temporarily turn off the S2 engine
  sf::sf_use_s2(FALSE) %>% suppressMessages()
  
  # Make point sf
  pt <- st_sf(
    geometry = st_sfc(st_point(c(lon, lat)), crs = 4326)
  )
  
  # Ensure both are in the same CRS
  message("Transforming data if needed ...")
  tictoc::tic("Completed.")
  if (st_crs(polygons) != st_crs(pt)) {
    pt <- st_transform(pt, st_crs(polygons))
  }
  tictoc::toc()
  
  # Spatial join: find the polygon that contains the point
  message("Performing join ...")
  tictoc::tic("Join complete.")
  result <- st_join(pt, polygons, join = st_intersects) %>% suppressMessages()
  tictoc::toc()
  
  # Drop geometry to just return attributes
  result_df <- result %>% st_drop_geometry()
  
  # turn S2 engine back on when done
  sf::sf_use_s2(TRUE) %>% suppressMessages()
  
  return(result_df)
}
