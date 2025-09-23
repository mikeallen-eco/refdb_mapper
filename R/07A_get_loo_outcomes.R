# compile leave-one-out data and apply taxonomy assignment rubrics & accuracy metrics

get_loo_outcomes <- function(marker_directories = dirname(refdb_cur_paths),
                             markers = markers,
                             refdb_nnd = NND_per_sp_within_refdb_seqs){
  
  # for(i in 1:length(marker_directories)){
  loo_list <- lapply(1:length(markers), function(i) {
    message("Processing marker: ", markers[i])
    
    # i=1 # TEMPORARY!
  loso_GB_compiled <- 
    compile_loo_ghostblaster_results(loo_path = paste0(marker_directories[i], "/loso/"))
  lospo_GB_compiled <- 
    compile_loo_ghostblaster_results(loo_path = paste0(marker_directories[i], "/lospo/"))
  
  refdb_nnd_df <- refdb_nnd %>%
    filter(marker %in% markers[i])
  
  loso_gb_outcomes <- loo_GB_outcomes(loso_GB_compiled) %>%
    mutate(marker = markers[i]) %>%
    left_join(refdb_nnd_df %>% mutate(assigned_mol_name = ununderscore(mol_name)) %>%
                select(assigned_mol_name, marker, nnd), 
              by = join_by(assigned_mol_name, marker))
  lospo_gb_outcomes <- loo_GB_outcomes(lospo_GB_compiled) %>%
    mutate(marker = markers[i])
  
  return(list(loso = loso_gb_outcomes,
              lospo = lospo_gb_outcomes))
  })
  
  names(loo_list) <- markers
  
}