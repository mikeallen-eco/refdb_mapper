

# Load packages
library(sf)
library(rnaturalearth)
library(dplyr)

# select a subset of US hydrobasins to analyze

subset_hydrobasins_US <- function(hydrobasin_map) {
  
# 1. Read your hydrobasins geopackage
hydrobasins <- st_read(hydrobasin_map)
  
# 2. Download US states from rnaturalearth
us_states <- ne_states(country = "United States of America", returnclass = "sf")

# 3. Keep only the 48 contiguous states
contig_us_states <- us_states %>%
  filter(!name %in% c("Alaska", "Hawaii", "Puerto Rico"))

# 4. Combine into a single US contiguous polygon
contig_us_boundary <- st_union(contig_us_states)

# 5. Reproject the boundary to match your hydrobasins CRS
#    (Hydrobasins are often in EPSG:4326, but confirm!)
hydrobasins_crs <- st_crs(hydrobasins)
contig_us_boundary <- st_transform(contig_us_boundary, hydrobasins_crs)

# 6. Identify which hydrobasins intersect the US boundary
#    This will give logical vector of TRUE/FALSE for each polygon
sf_use_s2(FALSE)
hydrobasins$in_US <- st_intersects(hydrobasins, contig_us_boundary, sparse = FALSE)[,1]

# 7. Optionally, filter to only basins inside US
hydrobasins_us <- hydrobasins %>%
  filter(in_US)

return(hydrobasins_us)

}