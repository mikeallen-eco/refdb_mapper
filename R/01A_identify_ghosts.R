
identify_ghosts <- function(hydrobasin_species = "~/Documents/mikedata/refdb_geo/hybas_L6_mammal_intersections_harmonized.csv", 
                            geotax_path = "data/geotax.csv", 
                            LOO_refdb_path){
  
# read in data
geotax <- read.csv(geotax_path) %>%
  dplyr::rename(sciname = orig_name)
h_data <- fread(hydrobasin_species)
LOO_refdb <- readDNAStringSet(LOO_refdb_path)

# get df of reference database
rn <- RDP_to_dataframe(LOO_refdb) %>%
  group_by(s) %>%
  mutate(n = length(s)) %>%
  ungroup() %>%
  mutate(join_name = gsub("_", " ", s))
refsp <- rn$join_name

# which hydrobasin species are ghosts
h_data_g <- h_data %>%
  left_join(select(geotax, sciname, gbif_name, ncbi_name), 
            by = join_by(sciname)) %>%
  mutate(ghost = case_when(!(sciname %in% refsp |
                             gbif_name %in% refsp |
                             ncbi_name %in% refsp) ~ 1,
                           TRUE ~ 0),
         join_name = case_when(ncbi_name %in% refsp ~ ncbi_name,
                               !(ncbi_name %in% refsp) & (gbif_name %in% refsp) ~ gbif_name,
                               !(ncbi_name %in% refsp) & !(gbif_name %in% refsp) ~ sciname)) %>%
  left_join(distinct(select(rn, join_name, n_seqs = n)), by = join_by(join_name)) %>%
  tidyr::replace_na(list(n_seqs = 0))

return(h_data_g)

}




