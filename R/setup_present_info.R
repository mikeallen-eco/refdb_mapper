library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(purrr, quietly = T, warn.conflicts = F)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "info|help")]
lapply(the_files_to_load, FUN = source)

markers = c("RiazVert1_12S", "MiMammalU_12S", "Vences_16S", "Mamm01_12S", "Taylor_16S")

final_sf <- readRDS("~/Documents/mikedata/refdb_mapper/final_hybas_data_sf_20251003.rds")
