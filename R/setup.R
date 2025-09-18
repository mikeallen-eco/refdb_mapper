library(data.table, quietly = T, warn.conflicts = F)
library(blastg, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(purrr, quietly = T, warn.conflicts = F)
library(stringdist, quietly = T, warn.conflicts = F)
library(Biostrings, quietly = T, warn.conflicts = F)
library(sf, quietly = T, warn.conflicts = F)
library(patchwork, quietly = T, warn.conflicts = F)
library(ape, quietly = T, warn.conflicts = F)
the_files <- list.files("R", full.names = T)
the_files_to_load <- the_files[grepl(the_files, pattern = "00|01|02|03|04|05|10|11|help")]
lapply(the_files_to_load, FUN = source)
hydrobasin_map <- read_sf("~/Documents/mikedata/refdb_mapper/hybas_L6_with_mammal_genus_richness.gpkg")
hydrobasin_species <- "~/Documents/mikedata/refdb_mapper/hybas_L6_mammal_intersections_harmonized.csv"
mol_tax <- "/Users/mikea/Documents/mikedata/mol_names/data/all_taxa_combined_taxonomy.csv"
tree_names <- "data/phyltax.csv"
manual_tax_refdb <- "data/refdb_mammals_manual_notes.tsv"
manual_tax_phyl <- "data/phyl_mammals_manual_notes.tsv"
refdb_harmonized_path <- "data/refdb_mammals_harmonized.csv"

refdb_RiazVert1_12S <- "~/Documents/mikedata/refdb_mapper/mammals_V5_12S/refdb_V512S_mammalia_midori265_tax20250609.fasta"
refdb_MiMammalU_12S <- "~/Documents/mikedata/refdb_mapper/mammals_MiMammalU_12S/refdb_MiMammalU12S_mammalia_midori265_tax20250609.fasta"
refdb_Vences_16S <- "~/Documents/mikedata/refdb_mapper/mammals_Vences_16S/refdb_V16S_mammalia_midori265_tax20250609.fasta"
refdb_Mamm01_12S <- "~/Documents/mikedata/refdb_mapper/mammals_Mamm01_12S/refdb_Mamm01_12S_mammalia_midori265_tax20250609.fasta"
refdb_Taylor_16S <- "~/Documents/mikedata/refdb_mapper/mammals_Taylor_16S/refdb_Taylor_16S_mammalia_midori265_tax20250609.fasta"

ncbi_extinct <- c("Homo_heidelbergensis", "Acratocnus_ye",
             "Arctodus_simus", "Bison_priscus",
             "Bison_schoetensacki", "Equus_dalianensis",
             "Camelus_knoblochi", "Equus_lambei",
             "Equus_ovodovi", "Bootherium_bombifrons",
             "Camelus_knoblochi", "Equus_hydruntinus",
             "Haringtonhippus_francisci", "Hippidion_saldiasi",
             "Homotherium_latidens", "Lasiopodomys_anglicus",
             "Hippidion_principale", "Homotherium_serum",
             "Hypnomys_morpheus", "Mammut_americanum",
             "Mammuthus_columbi", "Mammuthus_jeffersonii",
             "Mammuthus_primigenius", "Megaladapis_edwardsi",
             "Megaloceros_giganteus", "Megalonyx_jeffersonii",
             "Megatherium_americanum", "Miracinonyx_trumani",
             "Mylodon_darwinii", "Myotragus_balearicus", 
             "Nothrotheriops_shastensis", "Nesiotites_hidalgo", 
             "Panthera_spelaea",
             "Ursus_deningeri", "Ursus_ingressus", 
             "Ursus_kanivetz", "Ursus_kudarensis", "Ursus_spelaeus")

phyl_extinct <- c("Archaeolemur_majori", "Coelodonta_antiquitatis",
                  "Elephas_antiquus", "Elephas_cypriotes",
                  "Equus_capensis", "Eudorcas_rufina", 
                  "Homo_denisova", "Homo_neanderthalensis")
