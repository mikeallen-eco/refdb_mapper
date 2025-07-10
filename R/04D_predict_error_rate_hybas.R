predict_error_rate_hybas <- function(hybas_nnd_n_seqs,
                                     preds){
  
  pred_df <- lapply(hybas_nnd_n_seqs, function(x){
    
    preds_thresh99_i <- predict(
      preds$mod_thresh99_i,
      newdata = x,
      type = "response",
      re.form = NA
    )
    
    df <- x %>%
      mutate(thresh99_i = 100*preds_thresh99_i)
    
    return(df)
    
  }) %>%
    do.call(bind_rows, .)
  
  return(pred_df)
  
}
