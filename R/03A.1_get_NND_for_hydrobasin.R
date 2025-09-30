# get NND per species for a single hydrobasin and marker

get_NND_for_hydrobasin <- function(hydrobasin_refdb_info, 
                                   dist_full, 
                                   mol_phyl, 
                                   marker = "RiazVert1_12S", 
                                   HYBAS_ID_num = -1529923362,
                                   tree_path = "data/phyl.tre",
                                   extinct = c(ncbi_extinct, phyl_extinct)) {
  
  # read in the hydrobasin reference sequence count info & subset by marker
  hybas_sp_list <- hydrobasin_refdb_info %>%
    filter(HYBAS_ID %in% HYBAS_ID_num) %>%
    select(mol_name, n_seqs = marker)
  
  # subset to species with some sequences or zero sequences
  hybas_sp_list_no_zeros <- hybas_sp_list %>% filter(n_seqs > 0)
  hybas_sp_list_zeros    <- hybas_sp_list %>% filter(n_seqs == 0)
  
  # create lookup between MOL -> phylogeny
  phyl_map <- mol_phyl %>% filter(mol_name %in% hybas_sp_list$mol_name)
  
  # create empty list for results
  results <- list()
  
  # non-zeros: just compute within-group NND
  if (nrow(hybas_sp_list_no_zeros) > 1) {
    nonzeros <- phyl_map$phyl_name[match(hybas_sp_list_no_zeros$mol_name, phyl_map$mol_name)]
    nonzeros_no_na    <- nonzeros[!is.na(nonzeros)]
    nonzeros_mol    <- phyl_map$mol_name[match(hybas_sp_list_no_zeros$mol_name, phyl_map$mol_name)]
    nonzeros_mol_no_na    <- nonzeros_mol[!is.na(nonzeros_mol)]
    d_nonzeros <- apply(dist_full[nonzeros_no_na, nonzeros_no_na, drop = FALSE], 1, function(x) sort(x)[2])
    results$nonzeros <- tibble(mol_name = nonzeros_mol_no_na,
                               phyl_name = nonzeros_no_na,
                               nnd = d_nonzeros)
  }
  
  # zeros: min distance to nonzeros
  if (nrow(hybas_sp_list_zeros) > 0 && nrow(hybas_sp_list_no_zeros) > 0) {
    zeros    <- phyl_map$phyl_name[match(hybas_sp_list_zeros$mol_name, phyl_map$mol_name)]
    zeros_no_na    <- zeros[!is.na(zeros)]
    zeros_mol    <- phyl_map$mol_name[match(hybas_sp_list_zeros$mol_name, phyl_map$mol_name)]
    zeros_mol_no_na    <- zeros_mol[!is.na(zeros_mol)]
    nonzeros <- phyl_map$phyl_name[match(hybas_sp_list_no_zeros$mol_name, phyl_map$mol_name)]
    nonzeros_no_na    <- nonzeros[!is.na(nonzeros)]
    d_zeros  <- apply(dist_full[zeros_no_na, nonzeros_no_na, drop = FALSE], 1, min, na.rm = TRUE)
    results$zeros <- tibble(mol_name = zeros_mol_no_na,
                            phyl_name = zeros_no_na,
                            nnd = d_zeros)
  }
  
  # combine and attach HYBAS + marker info
  if(!isEmpty(results)){
  df <- bind_rows(results) %>%
    right_join(hybas_sp_list, by = "mol_name") %>%
    mutate(HYBAS_ID = HYBAS_ID_num)
  }else{df <- hybas_sp_list %>%
    mutate(HYBAS_ID = HYBAS_ID_num,
           nnd = NA)
  if(!"n_seqs" %in% colnames(df)){df$n_seqs <- NA}
  }
  
  # rename columns to include marker names
  colnames(df)[grepl("nnd|n_seqs", colnames(df))] <- paste0(marker, "_", colnames(df)[grepl("nnd|n_seqs", colnames(df))])
  
  df
}