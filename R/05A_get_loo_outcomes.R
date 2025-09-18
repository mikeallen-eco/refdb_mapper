# compile leave-one-out data and apply taxonomy assignment rubrics & accuracy metrics

get_loo_outcomes <- function(loso_gb_path = paste0(out_path, "loso/"),
                             lopso_gb_path = paste0(out_path, "lospo/"),
                             refdb_nnd = NND_per_sp_within_refdb){
  
  loso_GB_compiled <- compile_loo_ghostblaster_results(loo_path = loso_gb_path,
                                                       loo_refdb_nnd_df = refdb_nnd)
  
  lospo_GB_compiled <- compile_loo_ghostblaster_results(loo_path = lopso_gb_path,
                                                        loo_refdb_nnd_df = refdb_nnd)
  
  loso_gb_outcomes <- loo_GB_outcomes(loso_GB_compiled)
  lospo_gb_outcomes <- loo_GB_outcomes(lospo_GB_compiled)
  
  return(list(loso_gb_outcomes = loso_gb_outcomes,
              lospo_gb_outcomes = lospo_gb_outcomes))
  
}