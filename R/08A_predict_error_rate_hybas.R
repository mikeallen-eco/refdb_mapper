predict_error_rate_hybas <- function(hybas_info = hybas_nnd,
                                     preds = fits_blast97,
                                     markers = markers[1:5],
                                     assign_rubric){
  
  hybas_preds <- lapply(1:length(markers), function(i){

  message("Predicting error rates globally for marker: ", markers[i])    
  
  pred_df <- hybas_info %>%
    select(HYBAS_ID:phyl_name, contains(markers[i]))
  colnames(pred_df)[grepl(colnames(pred_df), pattern = "nnd|n_seqs")] <- c("nnd", "n_seqs")
  
  names(preds[[i]]) <- c("loso", "lospo") # to standardize in cases where loso_rdp etc. was used
  
  p_predictions = data.frame(
  loso_preds_i = predict(
    preds[[i]]$loso$i$mod,
    newdata = pred_df,
    type = "response",
    re.form = NA
  ),
  
  loso_preds_a = predict(
    preds[[i]]$loso$a$mod,
    newdata = pred_df,
    type = "response",
    re.form = NA
  ),
  
  loso_preds_c = predict(
    preds[[i]]$loso$c$mod,
    newdata = pred_df,
    type = "response",
    re.form = NA
  ),
  
  lospo_preds_i = predict(
    preds[[i]]$lospo$i$mod,
    newdata = pred_df,
    type = "response",
    re.form = NA
  ),
  
  lospo_preds_a = predict(
    preds[[i]]$lospo$a$mod,
    newdata = pred_df,
    type = "response",
    re.form = NA
  )
  )
  
  pred_df_final <- pred_df %>%
    bind_cols(p_predictions) %>%
    mutate(p_i = case_when(n_seqs %in% 0 ~ lospo_preds_i,
                           TRUE ~ loso_preds_i),
           p_a = case_when(n_seqs %in% 0 ~ lospo_preds_a,
                           TRUE ~ loso_preds_a),
           p_c = case_when(n_seqs %in% 0 ~ 0,
                           TRUE ~ loso_preds_c)) %>%
    select(-c(starts_with("loso"), starts_with("lospo")))
  
  colnames(pred_df_final)[grepl(colnames(pred_df_final), 
                                pattern = "nnd|n_seqs|p_")] <- paste0(markers[i], "_", 
                                                                      assign_rubric, "_",
                                                                      c("nnd", "n_seqs", "p_i",
                                                                        "p_a", "p_c"))

    return(pred_df_final)
  })
  
  names(hybas_preds) <- markers
  
  return(hybas_preds)

  }
