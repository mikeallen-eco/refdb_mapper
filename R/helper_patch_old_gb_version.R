patch_old_gb_version <- function(df) {
  
  # for testing
  # test_old1 <- 
  #   read.csv("~/Documents/mikedata/refdb_mapper/mammals_Vences_16S/loso/gd_1_Lepus_americanus_PQ049662.csv")
  # test_new1 <- 
  #   read.csv("~/Documents/mikedata/refdb_mapper/mammals_Mamm01_12S/loso/gd_1_Lepus_americanus_PQ049662.csv")
  # names(test_old1)
  # names(test_new1)[c(1:15, 21:31, 16:19, 32:37)]
  
  df_in <- df
  if("max_pident_loc_ecto" %in% names(df_in)){
  colnames(df_in) <- c("qseqid", "sseqid", "seq_species", "phyl_name", "group", "maxpident_blastg_local", 
                    "maxpident_blastg_all", "maxpident_blast", "cophenetic_from_blastg_local", 
                    "cophenetic_from_blastg_all", "cophenetic_from_blast_local", 
                    "cophenetic_from_blast_all", "in_refdb", "local", "in_phyl", 
                    "length", "mismatch", "gapopen", "qstart", "qend", "sstart", 
                    "send", "evalue", "bitscore", "skip_nohits", "badseq_index", 
                    "pred_max_pident_L", "pred_max_pident_M", "pred_max_pident_H", 
                    "blast_pred_metric", "BLAST_args", "min_length", "min_match", 
                    "max_treedepth", "n_seqs", "ref_seqs")
  }
 
  return(df_in) 
}