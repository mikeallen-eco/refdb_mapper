library(data.table)
library(blastg)
library(dplyr)
library(tidyr)
library(purrr)
library(stringdist)
library(Biostrings)
library(sf)
library(patchwork)
library(ape)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00|01A|02|03|04|05|10|11|help")]
lapply(the_files_to_load, FUN = source)
hydrobasin_map <- read_sf("~/Documents/mikedata/refdb_geo/hybas_L6_with_mammal_genus_richness.gpkg")
hydrobasin_species <- "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv"
mol_tax <- "/Users/mikea/Documents/mikedata/mol_names/data/all_taxa_combined_taxonomy.csv"
tree_names <- "data/phyltax.csv"
manual_tax_refdb <- "data/refdb_mammals_manual_notes.tsv"
refdb_harmonized_path <- "data/refdb_mammals_harmonized.csv"

extinct <- c("Homo_heidelbergensis", "Acratocnus_ye",
             "Arctodus_simus", "Bison_priscus",
             "Bison_schoetensacki", "Equus_dalianensis",
             "Camelus_knoblochi", "Equus_lambei",
             "Equus_ovodovi", "Bootherium_bombifrons",
             "Camelus_knoblochi", "Equus_hydruntinus",
             "Haringtonhippus_francisci", "Hippidion_saldiasi",
             "Homotherium_latidens", "Lasiopodomys_anglicus",
             "Hippidion_principale", "Homotherium_serum",
             "Hypnomys_morpheus")
