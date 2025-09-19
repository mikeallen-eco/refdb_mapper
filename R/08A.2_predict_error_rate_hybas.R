predict_error_rate_hybas <- function(hybas_nnd_n_seqs,
                                     preds){
  
  pred_df <- hybas_nnd_n_seqs %>%
    do.call(bind_rows, .)
    
    preds_i <- predict(
      preds$i$mod,
      newdata = pred_df,
      type = "response",
      re.form = NA
    )
    
    preds_a <- predict(
      preds$a$mod,
      newdata = pred_df,
      type = "response",
      re.form = NA
    )
    
    if("c" %in% names(preds)){
      preds_c <- predict(
        preds$c$mod,
        newdata = pred_df,
        type = "response",
        re.form = NA
      )
    }
    
    # predictions for future version of reference db w/ ≥ n seqs
    preds10_i <- predict(
      preds$i$mod,
      newdata = data.frame(nnd = pred_df$nnd, 
                           n_seqs = ifelse(pred_df$n_seqs < 10, 10, pred_df$n_seqs)),
      type = "response",
      re.form = NA
    )
    
    preds10_a <- predict(
      preds$a$mod,
      newdata = data.frame(nnd = pred_df$nnd, 
                           n_seqs = ifelse(pred_df$n_seqs < 10, 10, pred_df$n_seqs)),
      type = "response",
      re.form = NA
    )
    
    if("c" %in% names(preds)){
    preds10_c <- predict(
      preds$c$mod,
      newdata = data.frame(nnd = pred_df$nnd, 
                           n_seqs = ifelse(pred_df$n_seqs < 10, 10, pred_df$n_seqs)),
      type = "response",
      re.form = NA
    )
    }
    
    df <- pred_df %>%
      mutate(preds_i = 100*preds_i,
             preds10_i = 100*preds10_i,
             preds_a = 100*preds_a,
             preds10_a = 100*preds10_a)
    
    if("c" %in% names(preds)){
    df <- df %>%
      mutate(preds_c = 100*preds_c,
             preds10_c = 100*preds10_c)
    }
    
    return(df)
  
}
