# UNFINISHED!

# get evolutionary (cophenetic) distance for each species in a list of species (MOL)

get_NND_per_sp_within_list <- function(
             species_list, # vector of species names to loop through (from a MOL list)
             tree = "data/phylogeny_broader_species_pool.tre",
             phyltax = "data/phyltax.csv",
             broader_splist_tax_path = "data/geotax.csv",
             verbose = FALSE) {
  
  # load libraries and functions
  library(Biostrings, warn.conflicts = FALSE)
  library(dplyr)
  library(tidyr)
  library(ape)
  
  # read in the phylogenetic tree
  phyl_tree <- read.tree(tree) # combined tree covering all vertebrates

  dist.mat1 <- cophenetic.phylo(phyl_tree)
  row.names(dist.mat1) <- gsub("\\.", "-", row.names(dist.mat1)) # fix species with hyphens that got changed to periods
  colnames(dist.mat1) <- gsub("\\.", "-", colnames(dist.mat1)) # fix species with hyphens that got changed to periods
  
  # read in the broader SDM-based regional species list and the local SDM-based list
  broader_splist_tax <- read.csv(broader_splist_tax_path) %>%
    filter(!grepl(orig_name, pattern = " spp"))
  local_splist_tax <- broader_splist_tax %>%
    filter(orig_name %in% local_MOL_sp_list)
  
  ## get phylogeny names with taxonomy, subset to broader region
  # TODO: clean up this section if possible (maybe not necessary to create all these lists here)
  broader_allnames <- sort(unique(c(broader_splist_tax$orig_name, 
                                    broader_splist_tax$ncbi_name, 
                                    broader_splist_tax$gbif_name)))
  broader_allnames <- broader_allnames[!broader_allnames %in% "NA"]
  
  local_allnames <- sort(unique(c(local_splist_tax$orig_name, 
                                  local_splist_tax$ncbi_name, 
                                  local_splist_tax$gbif_name)))
  local_allnames <- local_allnames[!local_allnames %in% "NA"]
  
  # read in the phylogeny taxonomy lookup data, output of align_tax script
  phyltax_df_orig <- read.csv(phyltax) %>%
    filter(orig_name %in% broader_allnames | ncbi_name %in% broader_allnames | gbif_name %in% broader_allnames) %>%
    mutate(seq_species = case_when(!is.na(ncbi_name) ~ ncbi_name,
                                   TRUE ~ gbif_name)) %>%
    dplyr::rename(phyl_name = orig_name) %>% # renaming orig_name to make it clear it is the name used in the phylogeny
    # TODO: C. rufus snuck in as GBIF considers it to be preferably called Canis lupus; fix in align_tax script
    filter(!phyl_name %in% "Canis rufus") %>%
    # TODO: fix blank class in align_tax script for certain fish; a problem with GBIF class names
    mutate(class = case_when(order %in% c("Siluriformes", 
                                          "Characiformes",
                                          "Gymnotiformes",
                                          "Perciformes",
                                          "Cyprinodontiformes",
                                          "Clupeiformes",
                                          "Pleuronectiformes",
                                          "Beloniformes",
                                          "Synbranchiformes") ~ "Actinopteri",
                             TRUE ~ class)) %>%
    # adding this to ensure all names included in phyltax.csv are in fact on the tree
    filter(phyl_name %in% gsub("_", " ", phyl_tree$tip.label))
  
  allsp_results_list <- list()
  for (i in 1:length(species_list)) {
  message("Getting evgap for ", i, " of ", length(species_list), ": ", species_list[i])
  
    # gather info on focal NCBI species name
    focal.sp <- data.frame(ncbi_name = gsub("_", " ", species_list[i])) %>%
      left_join(select(phyltax_df_orig, phyl_name, 
                       gbif_name, ncbi_name,
                       seq_species),
                by = join_by(ncbi_name)) %>%
      slice_head(n = 1) %>% # in case there are >1 phyl_name per ncbi_name
      mutate(phyl_name2 = gsub(" ", "_", phyl_name))
    
    # check if focal species is in phylogeny
    if(is.na(focal.sp$phyl_name) | !focal.sp[1,]$seq_species %in% phyltax_df_orig$seq_species){
      if(verbose %in% TRUE){message("NOT FOUND IN PHYLOGENY!")}
      final_df <- data.frame(target_sp = species_list[i],
                             closest_rel_in_list = NA,
                             evgap_in_list = NA,      
                             closest_rel_all = NA,
                             evgap_all = NA)
    }else{
      if(verbose %in% TRUE){message("Detected in phylogeny!")}
    
    dist.mat <- dist.mat1 %>%
      as.data.frame() %>%
      dplyr::select(focal.sp$phyl_name2) %>%
      dplyr::rename(dist = 1) %>%
      arrange(dist) %>%
      tibble::rownames_to_column() %>%
      dplyr::rename(phyl_name = 1) %>%
      mutate(phyl_name2 = phyl_name,
             phyl_name = gsub("_", " ", phyl_name)) %>%
      slice_head(n = 1000) %>%
      left_join(select(phyltax_df_orig, phyl_name, 
                       gbif_name, ncbi_name,
                       seq_species),
                by = join_by(phyl_name)) %>%
      # add columns indicated whether in local species list and reference database
      mutate(in_list = case_when(ncbi_name %in% gsub("_", " ", species_list) ~ 1, 
                           TRUE ~ 0)) 

      dist.mat.final.in.list <- dist.mat %>%
        filter(in_list %in% 1) %>%
        arrange(dist) %>%
        slice_head(n = 2)
      
      dist.mat.final.all <- dist.mat %>%
        arrange(dist) %>%
        slice_head(n = 2)
      
      final_df <- data.frame(target_sp = species_list[i],
                             closest_rel_in_list = dist.mat.final.in.list$seq_species[2],
                             evgap_in_list = dist.mat.final.in.list$dist[2],
                             closest_rel_all = dist.mat.final.all$seq_species[2],
                             evgap_all = dist.mat.final.all$dist[2])
      
    } # end if not found in phylogeny
    
    allsp_results_list[[i]] <- final_df
  } # end species loop

  allsp_results <- do.call(rbind, allsp_results_list)
  
  return(allsp_results)
  
} # end function
