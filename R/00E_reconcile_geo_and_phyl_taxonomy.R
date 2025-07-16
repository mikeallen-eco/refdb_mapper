
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

none <- gbif_to_gbif %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match) 

mol_mam <- read.csv("~/Documents/mikedata/GhostBLASTer/2025_06/geotax_creation/MOL_MammaliaTaxonomy_v2.3_LF.csv")

none_list <- list()

for(i in 1:nrow(none)){

df <- mol_mam %>%
    filter(Accepted %in% none$geo_name[i]) %>%
    filter(Synonym %in% phyl_tax$phyl_name) %>%
  select(geo_name = Accepted, phyl_name = Synonym) %>%
  distinct() %>%
  group_by(geo_name) %>%
  slice_head(n = 1) %>%
  ungroup()

if(nrow(df) %in% 0){df <- data.frame(geo_name = none$geo_name[i], phyl_name = NA)}

none_list[[i]] <- df
  
}

none_df <- none_list %>%
  do.call(bind_rows, .) %>%
  mutate(match = "geo_to_syn")

none2 <- none_df %>%
  filter(is.na(phyl_name)) %>%
  select(-phyl_name, -match)


         in_phyl = case_when(!is.na(phyl_name) ~ 1,
                             TRUE ~ 0)) %>%
  select(order, geo_name, ncbi_name, gbif_name, phyl_name, in_phyl)

}