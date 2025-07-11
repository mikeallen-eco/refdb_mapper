get_preds_loso_lospo <- function(assign = "BLAST",
                                 accuracy_metric = "thresh98",
                                 outcomes){
  
  preds_gb_loso <- fit_model(assign_method = assign, loo_method = "LOSO", 
                             acc_metric = accuracy_metric, loo_outcomes_df = outcomes$loso_gb_outcomes)
  
  preds_gb_lospo <- fit_model(assign_method = assign, loo_method = "LOSpO", 
                              acc_metric = accuracy_metric, loo_outcomes_df = outcomes$lospo_gb_outcomes)
  
  return(list(preds_gb_loso = preds_gb_loso,
              preds_gb_lospo = preds_gb_lospo))
  
}