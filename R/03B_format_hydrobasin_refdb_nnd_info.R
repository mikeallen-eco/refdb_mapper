# combine NND and ref db info per hydrobasin into one df

format_hydrobasin_refdb_nnd_info <- function(hydrobasin_refdb_info, hybas_nnd){
  hybas_nnd <- hybas_nnd %>%
    do.call(bind_rows, .) %>%
    select(HYBAS_ID, mol_name, phyl_name, nnd)
  
  hydrobasin_refdb_nnd_info <- hydrobasin_refdb_info %>%
    left_join(hybas_nnd,
              by = join_by(HYBAS_ID, mol_name)) %>%
    select(HYBAS_ID, order, family, mol_name, phyl_name, nnd, contains("12S"), contains("16S")) %>%
    arrange(HYBAS_ID, order, nnd)
  
  return(hydrobasin_refdb_nnd_info)
}
