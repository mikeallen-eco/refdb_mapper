# compile leave-one-out data and apply taxonomy assignment rubrics & accuracy metrics

get_loo_outcomes <- function(marker_directories = dirname(refdb_cur_paths),
                             markers = markers,
                             refdb_harmonized_path = refdb_harmonized_path,
                             refdb_nnd = error_model_predictor_data){
  
  # for(i in 1:length(marker_directories)){
  loo_list <- lapply(1:length(markers), function(i) {
    message("Processing marker: ", markers[i])
    
  loso_GB_compiled <- 
    compile_loo_ghostblaster_results(loo_path = paste0(marker_directories[i], "/loso/"),
                                     refdb_harmonized_path = refdb_harmonized_path)
  
  lospo_GB_compiled <- 
    compile_loo_ghostblaster_results(loo_path = paste0(marker_directories[i], "/lospo/"),
                                     refdb_harmonized_path = refdb_harmonized_path)
  
  refdb_nnd_df <- refdb_nnd %>% mutate(true_mol_name = ununderscore(mol_name)) %>%
    filter(marker %in% markers[i]) %>%
    select(true_mol_name, phyl_name, marker, nnd, n_seqs)
  
  loso_gb_outcomes <- loo_GB_outcomes(loso_GB_compiled) %>%
    mutate(marker = markers[i]) %>%
    left_join(refdb_nnd_df %>% select(true_mol_name, phyl_name, nnd, n_seqs), 
              by = join_by(true_mol_name))
  
  lospo_gb_outcomes <- loo_GB_outcomes(lospo_GB_compiled) %>%
    mutate(marker = markers[i]) %>%
    left_join(refdb_nnd_df %>% select(true_mol_name, phyl_name, nnd, n_seqs), 
              by = join_by(true_mol_name))
  
  return(list(loso = loso_gb_outcomes,
              lospo = lospo_gb_outcomes))
  })
  
  names(loo_list) <- markers
  
  return(loo_list)
  
}