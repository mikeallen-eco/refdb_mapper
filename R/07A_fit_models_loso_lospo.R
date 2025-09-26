fit_models_loso_lospo <- function(assign_rubric = "thresh98",
                                 markers = markers,
                                 outcomes = mdat){
  
  preds_loso <- fit_models_one_loo_type(loo_method = "LOSO", 
                                 assign_rubric = assign_rubric,
                                 markers = markers,
                                 loo_outcomes_df = outcomes)
  
  preds_lospo <- fit_models_one_loo_type(loo_method = "LOSpO", 
                                  assign_rubric = assign_rubric,
                                  markers = markers,
                                  loo_outcomes_df = outcomes)
  
  # reorganize error model list structure to match: reorg$marker1$loso$i, reorg$marker1$lospo$i
  reorg <- reorganize_preds(preds_loso, preds_lospo)
  
  return(reorg)
  
}
