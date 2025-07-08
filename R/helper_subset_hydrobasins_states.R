


# Load packages
library(sf)
library(rnaturalearth)
library(dplyr)

# select a subset of US hydrobasins to analyze

subset_hydrobasins_states <- function(hydrobasin_map,
                                      states = c(
                                        "New Jersey",
                                        "Delaware",
                                        "Pennsylvania",
                                        "New York",
                                        "Massachusetts",
                                        "Rhode Island",
                                        "Maryland",
                                        "Connecticut",
                                        "Vermont",
                                        "New Hampshire",
                                        "Maine"
                                      )) {
  # 1. Read your hydrobasins geopackage
  hydrobasins <- st_read(hydrobasin_map)
  
  # 2. Download US states from rnaturalearth
  us_states <- ne_states(country = "United States of America", returnclass = "sf")
  
  # 3. Keep only the selected states
  us_states <- us_states %>%
    filter(name %in% states)
  
  # 4. Combine into a single US contiguous polygon
  region_boundary <- st_union(us_states)
  
  # 5. Reproject the boundary to match your hydrobasins CRS
  #    (Hydrobasins are often in EPSG:4326, but confirm!)
  hydrobasins_crs <- st_crs(hydrobasins)
  region_boundary <- st_transform(region_boundary, hydrobasins_crs)
  
  # 6. Identify which hydrobasins intersect the US boundary
  #    This will give logical vector of TRUE/FALSE for each polygon
  sf_use_s2(FALSE)
  hydrobasins$in_states <- st_intersects(hydrobasins, regions_boundary, sparse = FALSE)[, 1]
  
  # 7. Optionally, filter to only basins inside US
  hydrobasins_region <- hydrobasins %>%
    filter(in_states)
  
  return(hydrobasins_region)
  
}