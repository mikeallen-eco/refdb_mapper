

# Load packages
library(sf)
library(rnaturalearth)
library(dplyr)

# select a subset of US hydrobasins to analyze

subset_hydrobasins_countries <- function(hydrobasin_map,
                                         countries = "United States of America",
                                         exclude_states = c("Alaska", "Hawaii", "Puerto Rico")) {
  
# 1. Read your hydrobasins geopackage
hydrobasins <- st_read(hydrobasin_map)
  
# 2. Download US states from rnaturalearth
states <- ne_states(country = countries, returnclass = "sf")

# 3. Exclude states as indicated by exclude_states argument
contig_states <- states %>%
  filter(!name %in% exclude_states)

# 4. Combine into a single US contiguous polygon
contig_boundary <- st_union(contig_states)

# 5. Reproject the boundary to match your hydrobasins CRS
#    (Hydrobasins are often in EPSG:4326, but confirm!)
hydrobasins_crs <- st_crs(hydrobasins)
contig_oundary <- st_transform(contig_boundary, hydrobasins_crs)

# 6. Identify which hydrobasins intersect the US boundary
#    This will give logical vector of TRUE/FALSE for each polygon
sf_use_s2(FALSE)
hydrobasins$in_US <- st_intersects(hydrobasins, contig_boundary, sparse = FALSE)[,1]

# 7. Optionally, filter to only basins inside US
hydrobasins_us <- hydrobasins %>%
  filter(in_US)

return(hydrobasins_us)

}