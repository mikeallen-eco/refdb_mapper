harmonize_with_backbone <- function(query = mol_format,
                                    backbone = phyl_synonyms,
                                    manual_tax = NULL,
                                    fuzzy_threshold = 0.97) {
  
  # --- Step 1 - exact matches between query full_sci_name and backbone Accepted

  message("Step 1 - exact matches between query full_sci_name and backbone Accepted ...")
  
  # exact scientific name match (including subspecies designation where it is included)
  # match_type 1: direct match between query$full_sci_name and backbone$Accepted
  exact_sci <- query %>%
    select(uid, full_sci_name) %>%
    left_join(
      backbone %>% select(full_sci_name = Accepted) %>%
        mutate(match_type = "1 - exact_fullsci_w_ssp") %>%
        distinct(),
      by = join_by(full_sci_name)
    ) %>%
    mutate(bb_name = full_sci_name) %>%
    filter(match_type == "1 - exact_fullsci_w_ssp")

  message("   Matches found: ", nrow(exact_sci), " of ", nrow(query), " total\n")
  
  
  # --- Step 2 - exact matches between query full_sci_name and backbone Synonym
  
  message("Step 2 - exact matches between query full_sci_name and backbone Synonym ...")
  
  # exact synonym scientific name match (including query subspecies designations where given)
  # match_type 2: direct match between query$full_sci_name and backbone$Synonym
           
    exact_syn <- query %>%
    select(uid, full_sci_name) %>%
    filter(!uid %in% exact_sci$uid) %>%
    left_join(
      backbone %>% select(full_sci_name = Synonym, bb_name = Accepted) %>%
        mutate(match_type = "2 - exact_syn_w_ssp") %>%
        distinct(),
      by = join_by(full_sci_name),
      relationship = "many-to-many"
    ) %>%
    filter(match_type == "2 - exact_syn_w_ssp") %>%
    mutate(dup = 1 * duplicated(uid))
  
  dup_uids <- exact_syn %>%
    filter(dup %in% 1) %>%
    select(uid) %>%
    distinct() %>%
    pull()
  
  exact_syn2 <- exact_syn %>%
    filter(!uid %in% dup_uids) %>%
    select(-dup)
  
  message("   Additional matches: ", nrow(exact_syn2), "\n")
  if(length(dup_uids) > 0){
  message(
    "   This excludes ",
    length(dup_uids),
    " uids matched to multiple backbone Accepted names via synonyms ..."
  )
  message("   ... these records will be matched manually.\n")
  }
  
  
  # --- Step 3 - exact matches between query genus_species and backbone Accepted (excluding ssp name)
  
  message("Step 3 - exact matches between query genus_species and backbone Accepted (excluding ssp name) ...")
  
  # exact scientific name match (regardless of query subspecies designation)
  # match_type 3: direct match between query$"Genus Species" and backbone$Accepted
  exact_sci_wo_ssp <- query %>%
    select(uid, genus_species, full_sci_name) %>%
    filter(!uid %in% exact_sci$uid, !uid %in% exact_syn$uid) %>%
    left_join(
      backbone %>% select(genus_species = Accepted, bb_name = Accepted) %>%
        mutate(match_type = "3 - exact_sci_wo_ssp") %>%
        distinct(),
      by = join_by(genus_species)
    ) %>%
    filter(match_type == "3 - exact_sci_wo_ssp") %>%
    mutate(dup = 1 * duplicated(uid)) %>%
    select(-genus_species, -dup)

  message("   Additional matches: ", nrow(exact_sci_wo_ssp), "\n")
  
  
  # --- Step 4 - exact matches between query genus_species and backbone Synonym (excluding ssp name)

  message("Step 4 - exact matches between query genus_species and backbone Synonym (excluding ssp name) ...")
  
  # exact synonym scientific name match (regardless of query subspecies designation)
  # match_type 4: direct match between query$genus_species and backbone$Synonym
  exact_syn_wo_ssp <- query %>%
    select(uid, genus_species, full_sci_name) %>%
    filter(!uid %in% exact_sci$uid,!uid %in% exact_syn$uid,!uid %in% exact_sci_wo_ssp$uid) %>%
    left_join(
      backbone %>% select(genus_species = Synonym, bb_name = Accepted) %>%
        mutate(match_type = "4 - exact_syn_wo_ssp") %>%
        distinct(),
      by = join_by(genus_species)
    ) %>%
    filter(match_type == "4 - exact_syn_wo_ssp") %>%
    select(-genus_species)
  
  message("   Additional matches: ", nrow(exact_syn_wo_ssp), "\n")
  
  final_after_exact <- bind_rows(exact_sci, exact_syn2, 
                         exact_sci_wo_ssp, exact_syn_wo_ssp)
  
  # --- Step 5 - FUZZY matches between query full_sci_name and backbone Accepted (> threshold)
  
  message("Step 5 - fuzzy matches between query full_sci_name and backbone Accepted ...")
  
  # FUZZY scientific name match (including subspecies designation where it is included)
  # match_type 5: fuzzy version of match_type 1 to correct misspellings
  
  # function to get best Jaro-Winkler match and its score
  best_jw_match <- function(query, candidates) {
    sims <- stringsim(query, candidates, method = "jw")
    best_idx <- which.max(sims)
    list(best_match = candidates[best_idx], score = sims[best_idx])
  }
  
  candidates <- backbone$Accepted  
  
  message("   Fuzzy matching with Jaro-Winkler score threshold ", 
          fuzzy_threshold, " ...")
  
  fuzzy1 <- query %>% 
    filter(!uid %in% final_after_exact$uid) %>%
    mutate(match_res = map(full_sci_name, ~ best_jw_match(.x, candidates)),
           full_sci_name_fuz = map_chr(match_res, "best_match"),
           jw_score = map_dbl(match_res, "score")) %>%
    select(-match_res)
  
  # applying step 1 again with the fuzzy matches
  fuzzy_sci <- fuzzy1 %>%
    mutate(full_sci_name_orig = full_sci_name) %>%
    mutate(full_sci_name = case_when(jw_score >= fuzzy_threshold ~ full_sci_name_fuz,
                                     TRUE ~ full_sci_name)) %>%
    select(uid, full_sci_name, full_sci_name_orig, jw_score) %>%
    left_join(
      backbone %>% select(full_sci_name = Accepted) %>%
        mutate(match_type = "5 - fuzzy_fullsci_w_ssp") %>%
        distinct(),
      by = join_by(full_sci_name)
    ) %>%
    mutate(bb_name = full_sci_name) %>%
    filter(match_type == "5 - fuzzy_fullsci_w_ssp") %>%
    select(-full_sci_name) %>%
    dplyr::rename(full_sci_name = full_sci_name_orig)
  
  message("   Additional matches: ", nrow(fuzzy_sci), ".\n")
  
  # --- Step 6 - FUZZY mmatches between query full_sci_name and backbone Synonym (> threshold)
  
  message("Step 6 - fuzzy matches between query full_sci_name and backbone Synonym ...")
  
  # FUZZY scientific name match (including subspecies designation where it is included)
  # match_type 6: fuzzy version of match_type 2 to correct misspellings
  
  candidates_syn <- backbone$Synonym  
  
  message("   Fuzzy matching with Jaro-Winkler score threshold ", 
          fuzzy_threshold, " ...")
  
  fuzzy2 <- query %>% 
    filter(!uid %in% final_after_exact$uid) %>%
    filter(!uid %in% fuzzy_sci$uid) %>%
    filter(!uid %in% dup_uids) %>%
    mutate(match_res = map(full_sci_name, ~ best_jw_match(.x, candidates_syn)),
           full_sci_name_fuz = map_chr(match_res, "best_match"),
           jw_score = map_dbl(match_res, "score")) %>%
    select(-match_res)
  
  # applying step 2 again with the fuzzy matches
  fuzzy_syn <- fuzzy2 %>%
    mutate(full_sci_name_orig = full_sci_name) %>%
    mutate(full_sci_name = case_when(jw_score >= fuzzy_threshold ~ full_sci_name_fuz,
                                     TRUE ~ full_sci_name)) %>%
    select(uid, full_sci_name, full_sci_name_orig, jw_score) %>%
    left_join(
      backbone %>% select(full_sci_name = Synonym, bb_name = Accepted) %>%
        mutate(match_type = "6 - fuzzy_syn_w_ssp") %>%
        distinct(),
      by = join_by(full_sci_name)
    ) %>%
    filter(match_type == "6 - fuzzy_syn_w_ssp") %>%
    mutate(dup = 1 * duplicated(uid)) %>%
    select(-full_sci_name) %>%
    dplyr::rename(full_sci_name = full_sci_name_orig)
  
  dup_uids2 <- fuzzy_syn %>%
    filter(dup %in% 1) %>%
    select(uid) %>%
    distinct() %>%
    pull()
  
  fuzzy_syn2 <- fuzzy_syn %>%
    filter(!uid %in% dup_uids2) %>%
    select(-dup)
  
  message("   Additional matches: ", nrow(fuzzy_syn2))
  
  if(length(dup_uids2) > 0){
    message(
      "   This excludes ",
      length(dup_uids2),
      " uids matched to multiple backbone Accepted names via synonyms ..."
    )
    message("   ... these records will be matched manually.\n")
  }

  
  # --- Step 7 - FUZZY matches between query genus_species and backbone Accepted (> threshold)

  message("Step 7 - fuzzy matches between query genus_species and backbone Accepted (excluding ssp name) ...")  
  
  # FUZZY scientific name match (regardless of query subspecies designations)
  # match_type 7: fuzzy version of match_type 3 to correct misspellings
  
  candidates <- backbone$Accepted  
  
  message("   Fuzzy matching with Jaro-Winkler score threshold ", 
          fuzzy_threshold, ".")
  
  fuzzy3 <- query %>% 
    filter(!uid %in% final_after_exact$uid) %>%
    filter(!uid %in% unique(c(fuzzy_sci$uid, fuzzy_syn2$uid, dup_uids, dup_uids2))) %>%
    mutate(match_res = map(genus_species, ~ best_jw_match(.x, candidates)),
           genus_species_fuz = map_chr(match_res, "best_match"),
           jw_score = map_dbl(match_res, "score")) %>%
    select(-match_res)
  
  fuzzy_sci_wo_ssp <- fuzzy3 %>%
    mutate(genus_species = case_when(jw_score >= fuzzy_threshold ~ genus_species_fuz,
                                     TRUE ~ genus_species)) %>%
    select(uid, genus_species, full_sci_name, jw_score) %>%
    left_join(
      backbone %>% select(genus_species = Accepted, bb_name = Accepted) %>%
        mutate(match_type = "7 - fuzzy_sci_wo_ssp") %>%
        distinct(),
      by = join_by(genus_species)
    ) %>%
    filter(match_type == "7 - fuzzy_sci_wo_ssp") %>%
    select(-genus_species)
  
  message("   Additional matches: ", nrow(fuzzy_sci_wo_ssp), "\n")  
  
  
  # --- Step 8 - apply manual matches (including overrides) from manual_tax file
  
  message("Step 8 - applying manual matches ...")
  
  # match_type 8 (or 8.1-8.7 for manual overrides of match_types 1-7)
  
  if (is.null(manual_tax)) {
    message("   No manual_tax file specified. Skipping manual step for now.\n")
    manual <- data.frame(uid = NA, bb_name = NA, bb_name2 = NA, bb_name3 = NA, bb_name4 = NA)
  } 
  
  if (!is.null(manual_tax)) {
    
    message("   Reading manual_tax file ...")
    
    if (file.exists(manual_tax)) {
      man <- data.table::fread(manual_tax, header = T, sep = "\t")
    message("   Found ", nrow(man %>% filter(!is.na(bb_name))), " non-blank manual matches ...")
    } else{
      stop("Path to manual_tax file is invalid.")
    }
    
    manual <- query %>%
      select(uid, full_sci_name) %>%
      left_join(man %>% select(uid, bb_name, bb_name2, bb_name3, bb_name4), by = join_by(uid)) %>%
      mutate(match_type = case_when(uid %in% exact_sci$uid ~ "8.1 - manual - overriding type 1 match",
                                    uid %in% exact_syn2$uid ~ "8.2 - manual - overriding type 2 match",
                                    uid %in% exact_sci_wo_ssp$uid ~ "8.3 - manual - overriding type 3 match",
                                    uid %in% exact_syn_wo_ssp$uid ~ "8.4 - manual - overriding type 4 match",
                                    uid %in% fuzzy_sci$uid ~ "8.5 - manual - overriding type 5 match",
                                    uid %in% fuzzy_syn2$uid ~ "8.6 - manual - overriding type 6 match",
                                    uid %in% fuzzy_sci_wo_ssp$uid ~ "8.7 - manual - overriding type 7 match",
                                    uid %in% c(dup_uids, dup_uids2) ~ "8.8 - manual - resolved case with >1 backbone Accepted name",
                                    TRUE ~ "8.0 - manual - no automatic matches"))
    
    # remove blank rows
    manual <- manual %>% filter(!is.na(bb_name),
                                !bb_name %in% "")
  }
  
  # --- Step 9 - Final compiling and reporting 
  
  # bind all dfs together, removing manual matches from any previous steps 
  # (i.e., override them with the manual matches)
  final <- exact_sci %>%
    bind_rows(exact_syn2, exact_sci_wo_ssp, exact_syn_wo_ssp) %>%
    bind_rows(fuzzy_sci, fuzzy_syn2, fuzzy_sci_wo_ssp) %>%
    filter(!uid %in% manual$uid) %>%
    bind_rows(manual) %>%
    filter(!is.na(uid)) %>%
    right_join(query, by = join_by(uid, full_sci_name)) %>%
    mutate(dup = 1 * duplicated(uid)) %>%
    mutate(match_type = case_when(is.na(bb_name) & uid %in% c(dup_uids, dup_uids2) ~
                                    "9.1 - still unmatched, (>1 Backbone Accepted)",
                                  is.na(bb_name) & !uid %in% c(dup_uids, dup_uids2) ~
                                    "9.2 - still unmatched (no automatic matches)",
                                  TRUE ~ match_type)) %>%
    arrange(match_type, uid) %>%
    mutate(fuzzy_score = jw_score) %>%
    select(-jw_score, -genus_species) %>%
    dplyr::rename(BB_Accepted = bb_name, BB_Accepted2 = bb_name2, 
                  BB_Accepted3 = bb_name3, BB_Accepted4 = bb_name4)
  
  if(sum(final$dup) > 0){warning("Duplicate UIDs detected.")}
  
  # remove duplicate column
  final <- final %>% select(-dup)
  
  message("   Here is the final breakdown of match types ...\n")
  
  print(table(final$match_type))
  
  message("\n   It is advisable to check automatic matches, especially types 3-7.")
  message("   There are ",
          sum(is.na(final$BB_Accepted)),
          " records still unmatched.\n"
  )
  message("   Edit the manual_tax file and iterate to manually assign and override the backbone Accepted names.")
    message(
      "   The manual_tax file should be tsv format with at least 3 columns:\n",
      "      - uid is the unique query id.\n",
      "      - bb_name is the Backbone Accepted name assigned.\n",
      "      - bb_name2 is a second Backbone Accepted name assigned if applicable (i.e., for splits).\n",
      "      - bb_name3 is a third Backbone Accepted name assigned if applicable (i.e., for splits).\n",
      "      - bb_name4 is a fourth Backbone Accepted name assigned if applicable (i.e., for splits).\n"
    )
  
    if (sum(is.na(final$BB_Accepted)) > 0) {
      message("\nHere are the first 20 records that are still unmatched...")
      print(
        final %>% 
          filter(is.na(BB_Accepted)) %>% 
          select(uid, full_sci_name) %>%
          slice_head(n = 20)
      )
    }
    
  if("g" %in% names(final)){final <- final %>% select(-g)}
  if("s" %in% names(final)){final <- final %>% select(-s)}
  if("ssp" %in% names(final)){final <- final %>% select(-ssp)}

  message("Done.")
  
  return(final)
  
}