predict_error_rate_hybas <- function(hybas_nnd_n_seqs,
                                     preds){
  
  pred_df <- hybas_nnd_n_seqs %>%
    do.call(bind_rows, .)
    
    preds_i <- predict(
      preds$mod,
      newdata = pred_df,
      type = "response",
      re.form = NA
    )
    
    # predictions for 
    preds10_i <- predict(
      preds$mod,
      newdata = data.frame(nnd = pred_df$nnd, 
                           n_seqs = ifelse(pred_df$n_seqs < 10, 10, pred_df$n_seqs)),
      type = "response",
      re.form = NA
    )
    
    df <- pred_df %>%
      mutate(preds = 100*preds_i,
             preds10 = 100*preds10_i)
    
    return(df)
  
}
