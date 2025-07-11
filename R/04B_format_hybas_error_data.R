format_hybas_error_data <- function(hybas_ndd_df = hybas_nnd,
                                    ref_path = paste0(out_path, "refdb_", db_name, ".fasta"),
                                    preds_list = preds_loso_lospo) {
  hybas_nnd_n_seqs <- add_n_seqs_to_hybas_ndd(hybas_nnd = hybas_ndd_df, refdb = ref_path)
  
  hybas_pred_loso_i <- predict_error_rate_hybas(hybas_nnd_n_seqs = hybas_nnd_n_seqs,
                                                preds = preds_list$preds_gb_loso) %>%
    dplyr::rename(preds_loso = preds)
  
  hybas_pred_lospo_i <- predict_error_rate_hybas(hybas_nnd_n_seqs = hybas_nnd_n_seqs,
                                                 preds = preds_list$preds_gb_lospo) %>%
    dplyr::rename(preds_lospo = preds) %>% select(-preds10)
  
  hybas_pred_i <- hybas_pred_loso_i %>%
    left_join(
      hybas_pred_lospo_i,
      by = join_by(
        order,
        geo_name,
        seq_species,
        ncbi_name,
        phyl_name,
        in_phyl,
        nnd,
        HYBAS_ID,
        n_seqs
      )
    ) %>%
    mutate(
      preds = case_when(n_seqs %in% 0 ~ preds_lospo, n_seqs > 0 ~ preds_loso),
      preds_ghosts = case_when(n_seqs %in% 0 ~ preds_lospo, TRUE ~ NA)
    )
  
  return(hybas_pred_i)
  
}