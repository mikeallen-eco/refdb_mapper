# create taxonomy lookup for MOL -> phylogeny name
mol_to_phyl <- function(phyl_harmonized_path = phyl_harmonized_path){ 
  
  mol_to_phyl1 <- read.csv(phyl_harmonized_path) %>%
    filter(!is.na(MOL_Accepted)) %>%
    mutate(mol_name = underscore(MOL_Accepted),
           phyl_name = underscore(full_sci_name)) %>%
    select(mol_name, phyl_name)
  
  mol_to_phyl2 <- read.csv(phyl_harmonized_path) %>%
    filter(!is.na(MOL_Accepted2)) %>%
    mutate(mol_name = underscore(MOL_Accepted2),
           phyl_name = underscore(full_sci_name)) %>%
    select(mol_name, phyl_name)
  
  mol_to_phyl3 <- read.csv(phyl_harmonized_path) %>%
    filter(!is.na(MOL_Accepted3)) %>%
    mutate(mol_name = underscore(MOL_Accepted3),
           phyl_name = underscore(full_sci_name)) %>%
    select(mol_name, phyl_name)
  
  mol_to_phyl <- bind_rows(mol_to_phyl1,
                           mol_to_phyl2,
                           mol_to_phyl3) %>%
    group_by(mol_name) %>%
    slice_head(n=1) %>% # take the top match where > 1 phylogeny tip exists per MOL name
    ungroup()
  
  return(mol_to_phyl)
}