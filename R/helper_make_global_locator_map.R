# make a global locator map from a set of coordinates
library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)

make_global_locator_map <- function(lon = -74.759831, lat = 40.568704) {
  
  # Make point sf
  pt <- st_sf(geometry = st_sfc(st_point(c(lon, lat)), crs = 4326))
  
  pt_coords <- st_coordinates(pt)
  
  # Define orthographic projection centered on Eastern USA
  ortho_proj <- paste0(
    "+proj=ortho +lat_0=",
    round(pt_coords[, 2]),
    " +lon_0=",
    round(pt_coords[, 1]),
    " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  )
  
  # Load world basemap from Natural Earth
  world <- ne_countries(scale = "medium", returnclass = "sf")
  world_valid <- st_make_valid(world)
  world_parts <- st_cast(world_valid, "POLYGON") %>% suppressWarnings()  # preserves attributes
  
  # Transform both the basemap and your spatial data
  world_parts_ortho <- st_transform(world_parts, crs = ortho_proj)
  pt_ortho <- st_transform(pt, crs = ortho_proj)
  
  # Create circular mask using a buffer around a center point
  circle_mask <- st_point(c(0, 0)) |>
    st_sfc(crs = ortho_proj) |>
    st_buffer(dist = 6.4e6)  # ~6400 km = Earth's radius in meters
  
  # Plot
  (
    plot <- ggplot() +
      geom_sf(
        data = circle_mask,
        fill = "gray90",
        color = NA
      ) +  # background globe color
      geom_sf(
        data = world_parts_ortho,
        fill = "gray50",
        color = "white",
        size = 0.2
      ) +
      geom_sf(
        data = pt_ortho,
        color = "red",
        alpha = 1,
        size = 9,
        shape = "*"
      ) +
      coord_sf(crs = ortho_proj, datum = NA) +
      theme_void()
  )
  
  return(plot)
}