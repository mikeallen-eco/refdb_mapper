# function to get synonyms from rgbif

library(rgbif)
library(dplyr)

build_rgbif_synonym_database <- function(species) {
  
  get_synonyms <- function(name) {
    name_in <- ununderscore(name)
    
    # Try to get GBIF backbone info
    bb <- tryCatch(name_backbone(name_in), error = function(e) NULL)
    
    if (is.null(bb) || is.null(bb$usageKey)) {
      # No match in GBIF
      return(data.frame(
        Accepted = underscore(name_in),
        Synonym  = underscore(name_in),
        stringsAsFactors = FALSE
      ))
    }
    
    # Accepted name: prefer genus+species if available
    acc_name <- if (!is.null(bb$species)) {
      paste(bb$species)
    } else {
      bb$scientificName
    }
    acc_name <- underscore(acc_name)
    
    # Get synonyms safely
    syns <- tryCatch(name_usage(bb$usageKey, data = "synonyms"),
                     error = function(e) NULL)
    
    if (!is.null(syns) && !is.null(syns$data) && nrow(syns$data) > 0) {
      syns_sp_uni <- underscore(unique(sapply(
        strsplit(syns$data$scientificName, " "), 
        function(x) paste(x[1:2], collapse = " ")
      )))
      
      syns_df <- data.frame(
        Accepted = underscore(name_in),
        Synonym  = unique(c(acc_name, syns_sp_uni)),
        stringsAsFactors = FALSE
      )
    } else {
      syns_df <- data.frame(
        Accepted = underscore(name_in),
        Synonym  = acc_name,
        stringsAsFactors = FALSE
      )
    }
    return(syns_df)
  }
  
  synonyms_list <- lapply(seq_along(species), function(i) {
    message(i, "/", length(species), ": ", species[i])
    get_synonyms(species[i])
  })
  names(synonyms_list) <- species
  
  bind_rows(synonyms_list)
}


# library(rgbif)
# 
# build_rgbif_synonym_database <- function(species) {
#   
#   get_synonyms <- function(name) {
#     name <- ununderscore(name)
#     
#     # Get GBIF backbone info
#     bb <- name_backbone(name)
#     acc_name <- if (!is.null(bb$species)) {
#       paste(bb$species)
#     } else {
#       bb$scientificName
#     }
#     
#     key <- bb$usageKey
#     syns <- name_usage(key, data = "synonyms")
#     if (!is.null(syns$data) & nrow(syns$data) > 0) {
#       syns_sp_uni <- underscore(unique(sapply(
#         strsplit(syns$data$scientificName, " "), 
#         function(x) paste(x[1:2], collapse = " ")
#       )))
#       syns_df <- data.frame(Accepted = underscore(name), 
#                             Synonym = underscore(unique(c(acc_name, syns_sp_uni))))
#     } else {
#       syns_df <- data.frame(Accepted = underscore(name), Synonym = underscore(acc_name))
#     }
#     return(syns_df)
#   }
# 
# synonyms_list <- lapply(species, function(x) {message(which(species == x));get_synonyms(x)})
# names(synonyms_list) <- species
# synonyms_list
# 
# do.call(bind_rows, synonyms_list)
# }