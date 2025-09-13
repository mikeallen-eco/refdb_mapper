format_hybas_error_data <- function(hybas_nnd_df = hybas_nnd,
                                    ref_path = paste0(out_path, "refdb_", db_name, ".fasta"),
                                    preds_list = preds_loso_lospo) {
  hybas_nnd_n_seqs <- add_n_seqs_to_hybas_nnd(hybas_nnd = hybas_nnd_df, refdb = ref_path)
  
  hybas_pred_loso <- predict_error_rate_hybas(hybas_nnd_n_seqs = hybas_nnd_n_seqs,
                                                preds = preds_list$preds_gb_loso) %>%
    dplyr::rename_with(~ sub("^preds", "loso", .x))
  
  hybas_pred_lospo <- predict_error_rate_hybas(hybas_nnd_n_seqs = hybas_nnd_n_seqs,
                                                 preds = preds_list$preds_gb_lospo) %>%
    dplyr::rename_with(~ sub("^preds", "lospo", .x)) %>%
    dplyr::select(-lospo10_i, -lospo10_a)
  
  hybas_pred <- hybas_pred_loso %>%
    left_join(
      hybas_pred_lospo,
      by = join_by(
        order,
        geo_name,
        seq_species,
        ncbi_name,
        gbif_name,
        phyl_name,
        # in_phyl,
        nnd,
        HYBAS_ID,
        n_seqs
      )
    ) %>%
    mutate(
      preds_i = case_when(n_seqs %in% 0 ~ lospo_i, n_seqs > 0 ~ loso_i),
      preds_a = case_when(n_seqs %in% 0 ~ lospo_a, n_seqs > 0 ~ loso_a),
      preds_c = case_when(n_seqs %in% 0 ~ 0, n_seqs > 0 ~ loso_c)) %>%
    dplyr::rename(preds10_i = loso10_i,
                  preds10_a = loso10_a,
                  preds10_c = loso10_c)
  
  return(hybas_pred)
  
}