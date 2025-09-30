# # get evolutionary nearest neighbor distance of each species within all hydrobasins in a gpkg
# 
# # load libraries and functions
# library(data.table)
# library(dplyr)
# library(tidyr)
# library(ape)
# 
# get_NDD_per_sp_all_hydrobasins <- function(hydrobasin_refdb_info = hydrobasin_refdb_info,
#                                            markers = markers,
#                                            tree_path = "data/phyl.tre",
#                                            mol_to_phyl_path = mol_to_phyl_harmonized_path,
#                                            extinct = c(ncbi_extinct, phyl_extinct)){
#   
#   # read in info on hydrobasins (excluding 16 without any species)
#   hydrobasins <- data.frame(HYBAS_ID = unique(hydrobasin_refdb_info$HYBAS_ID))
#   n_hybas <- length(hydrobasins$HYBAS_ID)
#   
#   hybas_nnd_list <- lapply(1:n_hybas, function(x){
# 
#     df_list <- list()
#     
#     for(i in seq_along(markers)){
#       message(markers[i], ": ", x, " of ", n_hybas)
#       
#     hybas_sp_list <- hydrobasin_refdb_info %>%
#       filter(HYBAS_ID %in% hydrobasins$HYBAS_ID[x]) %>%
#       select(mol_name, n_seqs = markers[i])
#     
#     hybas_sp_list_no_zeros <- hybas_sp_list %>%
#       filter(n_seqs > 0)
#     
#     hybas_sp_list_zeros <- hybas_sp_list %>%
#       filter(n_seqs %in% 0)
#     
#     if(nrow(hybas_sp_list_no_zeros) > 0){
#     
#     df_no_zeros <- get_NND_per_sp_within_list(sp_list = hybas_sp_list_no_zeros$mol_name,
#                                tree = tree_path,
#                                mol_to_phyl_harmonized_path = mol_to_phyl_path,
#                                extinct = extinct)
#     
#     # loop through each species with 0 seqs and get nearest neighbor distance for species that DO have sequences
#     df_zeros <- lapply(seq_along(hybas_sp_list_zeros$mol_name), function(j){
#       message("zeros list: ", j, " of ", length(hybas_sp_list_zeros$mol_name))
#       splist <- unique(c(hybas_sp_list_zeros$mol_name[j], hybas_sp_list_no_zeros$mol_name))
#       df_one_sp <- get_NND_per_sp_within_list(sp_list = splist,
#                                                  tree = tree_path,
#                                                  mol_to_phyl_harmonized_path = mol_to_phyl_path,
#                                                  extinct = extinct) %>%
#       filter(mol_name %in% hybas_sp_list_zeros$mol_name[j])
#     
#     return(df_one_sp)
#     }) %>% do.call(bind_rows, .)
#       
#     df <- df_no_zeros %>%
#       bind_rows(df_zeros) %>%
#       right_join(hybas_sp_list, by = join_by(mol_name)) %>%
#       mutate(HYBAS_ID = hydrobasins$HYBAS_ID[x]) %>%
#       select(HYBAS_ID, mol_name, phyl_name, nnd, n_seqs)
#     colnames(df)[grepl(colnames(df), pattern = "nnd|n_seqs")] <- paste0(markers[i], "_", colnames(df)[grepl(colnames(df), pattern = "nnd|n_seqs")])
#     
#     }else{
#       
#       mol_phyl <- read.csv(mol_to_phyl_path) %>% 
#         mutate(mol_name = underscore(full_sci_name),
#                phyl_name = underscore(BB_Accepted)) %>%
#         select(mol_name, phyl_name)
#       
#       df <- hybas_sp_list %>% 
#         left_join(mol_phyl, by = join_by(mol_name)) %>%
#         mutate(HYBAS_ID = hydrobasins$HYBAS_ID[x], nnd = NA) %>%
#         select(HYBAS_ID, mol_name, phyl_name, nnd, n_seqs)
#       colnames(df)[grepl(colnames(df), pattern = "nnd|n_seqs")] <- paste0(markers[i], "_", colnames(df)[grepl(colnames(df), pattern = "nnd|n_seqs")])
#     }
#     
#     df_list[[i]] <- df
#     }
#     
#     names(df_list) <- markers
#     
#     return(df_list)
#   })
#   
#   return(hybas_nnd_list)
#   
# }
