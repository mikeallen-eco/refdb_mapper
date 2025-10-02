library(purrr)

fit_models_loso_lospo <- function(assign_rubric = "thresh98",
                                 markers = markers,
                                 outcomes = mdat){
  
  preds_loso <- fit_models_one_loo_type(loo_method = "LOSO", 
                                 assign_rubric = assign_rubric,
                                 markers = markers,
                                 loo_outcomes_list = outcomes)
  
  preds_lospo <- fit_models_one_loo_type(loo_method = "LOSpO", 
                                  assign_rubric = assign_rubric,
                                  markers = markers,
                                  loo_outcomes_list = outcomes)
  
  # function to reorganize list to match: reorg$marker1$loso$i, reorg$marker1$lospo$i
  reorganize_preds <- function(preds_loso, preds_lospo) {
    # get union of marker names across both
    markers <- union(names(preds_loso), names(preds_lospo))
    
    # build reorganized list
    purrr::map(markers, function(m) {
      list(
        loso  = preds_loso[[m]],
        lospo = preds_lospo[[m]]
      )
    }) |> set_names(markers)
  }
  
  # reorganize list to match: reorg$marker1$loso$i, reorg$marker1$lospo$i
  reorg <- reorganize_preds(preds_loso, preds_lospo)
  
  return(reorg)
  
}
