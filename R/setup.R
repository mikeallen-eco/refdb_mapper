library(data.table)
library(ghostblaster)
library(dplyr)
library(Biostrings)
library(sf)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00_|00A|00B|00C|00D|01|02|03|04|05|06|help")]
lapply(the_files_to_load, FUN = source)
hydrobasin_map <- st_read("~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg")
