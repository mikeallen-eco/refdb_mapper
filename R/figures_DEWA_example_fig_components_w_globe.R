library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)

# Filter your feature data
dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

# Load world basemap from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define orthographic projection centered on Eastern USA
ortho_proj <- "+proj=ortho +lat_0=40 +lon_0=-75 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform both the basemap and your spatial data
world_ortho <- st_transform(world, crs = ortho_proj)
dewa_ortho <- st_transform(dewa_hybas_pred_map_sf, crs = ortho_proj)
hybas_ortho <- st_transform(hybas_pred_map_sf, crs = ortho_proj)

# Create circular mask using a buffer around a center point
circle_mask <- st_point(c(0, 0)) |>
  st_sfc(crs = ortho_proj) |>
  st_buffer(dist = 6.4e6)  # ~6400 km = Earth's radius in meters

# Plot
ggplot() +
  geom_sf(data = circle_mask, fill = "gray90", color = NA) +  # background globe color
  geom_sf(data = world_ortho, fill = "gray50", color = "white", size = 0.2) +
  geom_sf(data = hybas_ortho, fill = "gray50", color = "gray80", size = 0.2) +
  geom_sf(data = dewa_ortho, fill = "red", alpha = 0.7, color = NA) +
  coord_sf(crs = ortho_proj) +
  theme_void() +
  theme(text = element_text(size = 18))

# Save figure
ggsave("figures/conceptual/globe_map_wo_hybas.png", height = 6, width = 6, dpi = 400)
