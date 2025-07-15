library(data.table)
library(ghostblaster)
library(dplyr)
library(Biostrings)
library(sf)
library(ggplot2)
library(patchwork)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00|01|02|03|04|05|help")]
lapply(the_files_to_load, FUN = source)
out_path <- "~/Documents/mikedata/refdb_geo/mammals_Vences_16S/"
db_name <- "V16S_mammalia_midori265_tax20250609"