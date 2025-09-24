fit_models_loso_lospo <- function(assign_rubric = "thresh98",
                                 markers = markers[1:3],
                                 outcomes = mdat){
  
  preds_loso <- fit_models_one_loo_type(loo_method = "LOSO", 
                                 assign_rubric = assign_rubric,
                                 markers = markers,
                                 loo_outcomes_df = outcomes)
  
  preds_lospo <- fit_models_one_loo_type(loo_method = "LOSpO", 
                                  assign_rubric = assign_rubric,
                                  markers = markers,
                                  loo_outcomes_df = outcomes)
  
  return(list(loso = preds_loso,
              lospo = preds_lospo))
  
}
