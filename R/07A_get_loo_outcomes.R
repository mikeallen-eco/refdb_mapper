# compile leave-one-out data and apply taxonomy assignment rubrics & accuracy metrics

get_loo_outcomes <- function(marker_directories = dirname(refdb_cur_paths),
                             refdb_nnd = NND_per_sp_within_refdb){
  
  # for(i in 1:length(marker_directories)){
    i=1 # TEMPORARY!
  loso_GB_compiled <- compile_loo_ghostblaster_results(loo_path = paste0(marker_directories[i], "/loso/"),
                                                       loo_refdb_nnd_df = refdb_nnd)
  
  lospo_GB_compiled <- compile_loo_ghostblaster_results(loo_path = lopso_gb_path,
                                                        loo_refdb_nnd_df = refdb_nnd)
  
  loso_gb_outcomes <- loo_GB_outcomes(loso_GB_compiled)
  lospo_gb_outcomes <- loo_GB_outcomes(lospo_GB_compiled)
  
  return(list(loso_gb_outcomes = loso_gb_outcomes,
              lospo_gb_outcomes = lospo_gb_outcomes))
  
}