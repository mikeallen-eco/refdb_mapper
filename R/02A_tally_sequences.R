library(data.table)

tally_sequences <- function(hydrobasin_species = hydrobasin_species_path, 
                            refdb_cur_paths = refdb_cur_paths,
                            refdb_harmonized_path = refdb_harmonized_path,
                            extinct = ncbi_extinct){
  
# read in data

h_data <- fread(hydrobasin_species)
refdb_harmonized <- read.csv(refdb_harmonized_path)

# readDNAStringSet(refdb_cur_path)
refdb_cur <- lapply(refdb_cur_paths, function(x) Biostrings::readDNAStringSet(x) )

h_data_g_list <- lapply(refdb_cur, function(x) {
       
# get count of sequences for each MOL name
mol_ref1 <- RDP_to_dataframe(x) %>%
  filter(!s %in% extinct,
         !grepl(s, pattern = "_x_")) %>%
  group_by(s) %>%
  summarize(n = length(s), 
            .groups = "drop") %>%
  left_join(refdb_harmonized %>% select(s = uid, BB_Accepted),
            by = join_by(s))

mol_ref2 <- RDP_to_dataframe(x) %>%
  filter(!s %in% extinct,
         !grepl(s, pattern = "_x_")) %>%
  group_by(s) %>%
  summarize(n = length(s), 
            .groups = "drop") %>%
  left_join(refdb_harmonized %>% select(s = uid, BB_Accepted = BB_Accepted2),
            by = join_by(s))

mol_ref <- bind_rows(mol_ref1, mol_ref2) %>%
  filter(!is.na(BB_Accepted)) %>%
  group_by(BB_Accepted) %>%
  summarize(n = sum(n),
            .groups = "drop")

# which hydrobasin species are ghosts
h_data_g <- h_data %>%
  left_join(mol_ref %>% select(sciname = BB_Accepted, n_seqs = n), 
            by = join_by(sciname)) %>%
  tidyr::replace_na(list(n_seqs = 0))

return(h_data_g)

})

#  h_data_g_list is list of data.tables
h_merged <- Reduce(function(x, y) merge(x, y, 
                                        by = c("HYBAS_ID", "sciname", "order", "family", "genus", "species"),
                                        all = TRUE),
                   lapply(names(h_data_g_list), function(nm) {
                     dt <- copy(h_data_g_list[[nm]])
                     setnames(dt, old = c("n_seqs"), 
                              new = paste0(nm))
                     dt
                   })) %>%
  dplyr::rename(mol_name = sciname) %>%
  mutate(mol_name = underscore(mol_name)) %>%
  select(-genus, -species)

return(h_merged)

}