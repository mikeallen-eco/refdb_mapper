
reconcile_geo_to_phyl_taxonomy <- function(
    phyltax = "data/phyltax.csv",
    sp_list_tax = "data/geotax.csv",
    verbose = FALSE) {
  
# read in phylogenetic tree taxonomy
phyl_tax <- read.csv(phyltax) %>%
  filter(class %in% "Mammalia")

# read in input species list taxonomy
geo_tax <- read.csv(sp_list_tax) %>%
  dplyr::rename(geo_name = orig_name) %>%
  filter(class %in% "Mammalia") %>%
  mutate(ncbi_name = case_when(grepl(ncbi_name, pattern = "sp\\.| x ") ~ NA,
                               TRUE ~ ncbi_name))

geo_to_phyl <- geo_tax %>%
  left_join(phyl_tax %>% mutate(geo_name = phyl_name, match = "geo_to_phyl") %>% 
              select(geo_name, phyl_name, match) %>% distinct(), 
            by = join_by(geo_name))

geo_to_ncbi <- geo_to_phyl %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(geo_name = ncbi_name, match = "geo_to_ncbi") %>%
            select(geo_name, phyl_name, match) %>% distinct(), 
          by = join_by(geo_name))

geo_to_gbif <- geo_to_ncbi %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(geo_name = gbif_name, match = "geo_to_gbif") %>%
              select(geo_name, phyl_name, match) %>% distinct(), 
            by = join_by(geo_name))

ncbi_to_phyl <- geo_to_gbif %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(ncbi_name = phyl_name, match = "ncbi_to_phyl") %>%
              select(ncbi_name, phyl_name, match) %>% distinct(),
            by = join_by(ncbi_name))

ncbi_to_ncbi <- ncbi_to_phyl %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(match = "ncbi_to_ncbi") %>%
              filter(!is.na(ncbi_name)) %>%
              select(ncbi_name, phyl_name, match) %>% distinct(),
            by = join_by(ncbi_name))

ncbi_to_gbif <- ncbi_to_ncbi %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(ncbi_name = gbif_name, match = "ncbi_to_gbif") %>%
              filter(!is.na(ncbi_name)) %>%
              select(ncbi_name, phyl_name, match) %>% distinct(),
            by = join_by(ncbi_name))
    
gbif_to_phyl <- ncbi_to_gbif %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(gbif_name = phyl_name, match = "gbif_to_phyl") %>%
              filter(!is.na(gbif_name)) %>%
              select(gbif_name, phyl_name, match) %>% distinct(),
            by = join_by(gbif_name))

gbif_to_ncbi <- gbif_to_phyl %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(gbif_name = ncbi_name, match = "gbif_to_ncbi") %>%
              filter(!is.na(gbif_name)) %>%
              select(gbif_name, phyl_name, match) %>% distinct(),
            by = join_by(gbif_name))

gbif_to_gbif <- gbif_to_ncbi %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  left_join(phyl_tax %>% mutate(match = "gbif_to_gbif") %>%
              filter(!is.na(gbif_name)) %>%
              select(gbif_name, phyl_name, match) %>% distinct(),
            by = join_by(gbif_name))

syn_to_phyl1 <- gbif_to_gbif %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) 

mol_mam <- read.csv("~/Documents/mikedata/GhostBLASTer/2025_06/geotax_creation/MOL_MammaliaTaxonomy_v2.3_LF.csv")

syn_to_phyl_list <- list()

for(i in 1:nrow(syn_to_phyl1)){

df <- mol_mam %>%
    filter(Accepted %in% syn_to_phyl1$geo_name[i]) %>%
    filter(Synonym %in% phyl_tax$phyl_name) %>%
  select(geo_name = Accepted, phyl_name = Synonym) %>%
  distinct() %>%
  group_by(geo_name) %>%
  slice_head(n = 1) %>%
  ungroup()

if(nrow(df) %in% 0){df <- data.frame(geo_name = syn_to_phyl1$geo_name[i], phyl_name = NA)}

syn_to_phyl_list[[i]] <- df
  
}

syn_to_phyl <- syn_to_phyl_list %>%
  do.call(bind_rows, .) %>%
  mutate(match = "syn_to_phyl")

syn_to_ncbi1 <- syn_to_phyl %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match)

syn_to_ncbi_list <- list()

for(i in 1:nrow(syn_to_ncbi1)){
  
  df <- mol_mam %>%
    filter(Accepted %in% syn_to_ncbi1$geo_name[i]) %>%
    filter(Synonym %in% phyl_tax$ncbi_name) %>%
    select(geo_name = Accepted, phyl_ncbi_name = Synonym) %>%
    distinct() %>%
    group_by(geo_name) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if(nrow(df) %in% 0){df <- data.frame(geo_name = syn_to_ncbi1$geo_name[i], phyl_ncbi_name = NA)}
  
  syn_to_ncbi_list[[i]] <- df
  
}

syn_to_ncbi2 <- syn_to_ncbi_list %>%
  do.call(bind_rows, .) %>%
  mutate(match = "syn_to_ncbi")

phyl_ncbi_tmp <- phyl_tax %>%
  filter(!is.na(ncbi_name),
         ncbi_name %in% syn_to_ncbi2$phyl_ncbi_name) %>%
  select(phyl_name, phyl_ncbi_name = ncbi_name)

syn_to_ncbi <- syn_to_ncbi2 %>%
  left_join(phyl_ncbi_tmp, by = join_by(phyl_ncbi_name)) %>%
  select(geo_name, phyl_name, match)

# try to match geo synonyms with phyl gbif names
syn_to_gbif1 <- syn_to_ncbi %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match)

syn_to_gbif_list <- list()

for(i in 1:nrow(syn_to_gbif1)){
  
  df <- mol_mam %>%
    filter(Accepted %in% syn_to_gbif1$geo_name[i]) %>%
    filter(Synonym %in% phyl_tax$gbif_name) %>%
    select(geo_name = Accepted, phyl_gbif_name = Synonym) %>%
    distinct() %>%
    group_by(geo_name) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if(nrow(df) %in% 0){df <- data.frame(geo_name = syn_to_gbif1$geo_name[i], phyl_gbif_name = NA)}
  
  syn_to_gbif_list[[i]] <- df
  
}

syn_to_gbif2 <- syn_to_gbif_list %>%
  do.call(bind_rows, .) %>%
  mutate(match = "syn_to_gbif")

phyl_gbif_tmp <- phyl_tax %>%
  filter(!is.na(gbif_name),
         gbif_name %in% syn_to_gbif2$phyl_gbif_name) %>%
  select(phyl_name, phyl_gbif_name = gbif_name)

syn_to_gbif <- syn_to_gbif2 %>%
  left_join(phyl_gbif_tmp, by = join_by(phyl_gbif_name)) %>%
  select(geo_name, phyl_name, match)

geo_to_phyl_manual <- syn_to_gbif %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) %>%
  mutate(phyl_name = case_when(geo_name %in% c("Giraffa giraffa", 
                            "Giraffa tippelskirchi") ~ "Giraffa camelopardalis",
                            geo_name %in% "Sciurus meridionalis" ~ "Sciurus vulgaris",
                            geo_name %in% "Acomys chudeaui" ~ "Acomys cahirinus",
                            geo_name %in% "Balaenoptera ricei" ~ "Balaenoptera edeni",
                            TRUE ~ NA),
         match = "manual")
  
  no_matches <- geo_to_phyl_manual %>%
    filter(is.na(phyl_name)) %>%
    select(-match) %>%
    mutate(match = "no_match") %>%
    arrange(geo_name)
    
final_geo_to_phyl_df <- geo_to_phyl %>%
  bind_rows(geo_to_ncbi, geo_to_gbif, 
            ncbi_to_phyl, ncbi_to_ncbi, ncbi_to_gbif,
            gbif_to_phyl, gbif_to_ncbi, gbif_to_gbif,
            syn_to_phyl, syn_to_ncbi, syn_to_gbif,
            geo_to_phyl_manual) %>%
  filter(!is.na(phyl_name)) %>%
  bind_rows(no_matches) %>%
  group_by(geo_name) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::rename(phyl_match = match) %>%
  select(geo_name, phyl_name, phyl_match)

geo_tax_final <- geo_tax %>%
  left_join(final_geo_to_phyl_df, 
            by = join_by(geo_name))

write.csv(geo_tax_final, "data/geotax_phyl_mammals.csv", row.names = F)

}