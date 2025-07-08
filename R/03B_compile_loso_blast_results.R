# evaluate BLAST outcomes for leave-one-sequence-out data & add nearest cophenetic distance info

compile_loso_blast_results <- function(loso_path = paste0(out_path, "loso"),
                                       loso_refdb_ndd_df = loso_refdb_ndd) {
  
  loso <- read_csvs(loso_path, pattern_str = "csv")
  
  loso_list <- split(loso, loso$qseqid)
  
  loso_blast_df <- lapply(
    loso_list,
    FUN = function(x) {
      pick_top_match <- x %>%
        select(qseqid, seq_species, phyl_name, max_pident = max_pident_no_ecto) %>%
        filter(!is.na(max_pident)) %>%
        mutate(max_pident_c = as.character(max_pident)) %>%
        group_by(max_pident_c) %>%
        mutate(n = length(max_pident)) %>%
        slice_sample(n = 1) %>% # pick 1 at random in cases of identical max_pident
        ungroup() %>%
        arrange(desc(max_pident))
      
      df <- data.frame(
        qseqid = pick_top_match$qseqid[1],
        top_ncbi_name = pick_top_match$seq_species[1],
        next_ncbi_name = pick_top_match$seq_species[2],
        top_pident = sort(pick_top_match$max_pident, decreasing = T)[1],
        next_pident = sort(pick_top_match$max_pident, decreasing = T)[2],
        n_top = pick_top_match$n[1]
      ) %>%
        mutate(gap = top_pident - next_pident, tmp = qseqid) %>%
        separate(tmp, into = c("g", "s", "acc"), sep = "_") %>%
        mutate(true_ncbi_name = paste0(g, " ", s),
               gap = case_when(n_top > 1 ~ 0, TRUE ~ gap)) %>%
        select(-g, -s, acc) %>%
        select(qseqid,
               true_ncbi_name,
               top_ncbi_name,
               max_pident = top_pident,
               n_top,
               gap)
      
      return(df)
    }
  ) %>%
    do.call(bind_rows, .)
  
  loso_blast_df_filtered <- loso_blast_df %>%
    filter(true_ncbi_name %in% loso_refdb_ndd$ncbi_name) %>%
    left_join(loso_refdb_ndd %>% dplyr::rename(true_ncbi_name = ncbi_name),
              by = join_by(true_ncbi_name)) %>%
    select(-phyl_name) %>%
    mutate(blast97 = case_when((
      max_pident >= 97 &
        (max_pident - gap) < 97 &
        top_ncbi_name == true_ncbi_name
    ) ~ "correct",
    (
      max_pident >= 97 &
        max_pident - gap < 97 &
        top_ncbi_name != true_ncbi_name
    ) ~ "incorrect",
    !(max_pident >= 97 &
        (max_pident - gap) < 97) ~ "abstain"
    )) %>%
    mutate(blast98 = case_when((
      max_pident >= 98 &
        (max_pident - gap) < 98 &
        top_ncbi_name == true_ncbi_name
    ) ~ "correct",
    (
      max_pident >= 98 &
        (max_pident - gap) < 98 &
        top_ncbi_name != true_ncbi_name
    ) ~ "incorrect",
    !(max_pident >= 98 &
        (max_pident - gap) < 98) ~ "abstain"
    )) %>%
    mutate(blast99 = case_when((
      max_pident >= 99 &
        (max_pident - gap) < 99 &
        top_ncbi_name == true_ncbi_name
    ) ~ "correct",
    (
      max_pident >= 99 &
        (max_pident - gap) < 99 &
        top_ncbi_name != true_ncbi_name
    ) ~ "incorrect",
    !(max_pident >= 99 &
        (max_pident - gap) < 99) ~ "abstain"
    )) %>%
    mutate(ecotag = case_when(((100 - max_pident) < gap &
                                 top_ncbi_name == true_ncbi_name) ~ "correct",
                              ((100 - max_pident) < gap &
                                 top_ncbi_name != true_ncbi_name) ~ "incorrect",
                              !((100 -
                                   max_pident) < gap) ~ "abstain"
    ))
  
  return(loso_blast_df_filtered)
}