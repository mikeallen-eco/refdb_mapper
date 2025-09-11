# read in and format data
prepare_target_data_for_harmonization <- function(mol_tax = mol_tax,
                                                  mol_group = c("Mammals"),
                                                  tree_names = tree_names,
                                                  refdb_cur = refdb_cur,
                                                  extinct = extinct) {
  
  mol <- read.csv(mol_tax) %>%
    filter(taxa %in% mol_group)
  
  if (mol_group %in% "Mammals"){phyl_groups <- c("Mammalia")}else(phyl_groups <- NULL)
      
  phyl <- read.csv(tree_names) %>%
    select(uid = phyl_name2,
           full_sci_name = phyl_name,
           class) %>%
    mutate(tmp = full_sci_name) %>%
    separate(tmp, into = c("g", "s", "ssp"), sep = " ") %>%
    filter(!g %in% "X",
           class %in% phyl_groups) %>%
    mutate(genus_species = paste0(g, " ", s)) %>%
    suppressWarnings()
  
  refdb <- lapply(refdb_cur, function(x) Biostrings::readDNAStringSet(x) ) %>%
    do.call("c", .) %>%
    names() %>%
    as.data.frame() %>%
    dplyr::rename(header = 1) %>%
    separate(header, into = c("acc", "string"), sep = "\t") %>%
    separate(string, into = c("root", "k", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
    select(-acc, -root, -k) %>%
    distinct() %>%
    mutate(full_sci_name = gsub("_", " ", s)) %>%
    dplyr::rename(uid = s) %>%
    filter(!uid %in% extinct,
           !grepl(uid, pattern = "_x_")) %>%
    mutate(tmp = full_sci_name) %>%
    separate(tmp, into = c("g", "s", "ssp"), sep = " ") %>%
    mutate(genus_species = paste0(g, " ", s)) %>%
    suppressWarnings()
    
data_list <- list(mol = mol,
                  phyl = phyl,
                  refdb = refdb)
  
  return(data_list)
  
}
